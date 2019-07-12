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
C------------------------------------------------------------------------
C||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
C------------------------------------------------------------------------
C File: abacus/herqm3.F
C
C Written by JK+AO 2002
C
C-------------------------------------------------------------------------
C  /* Deck qm3inp */
      SUBROUTINE QM3INP(WORD)
C-------------------------------------------------------------------------
C
C INPUT FROM ABACUS TO QM3 CALCULATIONS 
C
#include "implicit.h"
#include "priunit.h"
#include "gnrinf.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "mxsymm.h"
#include "qm3.h"
#include "maxorb.h"
#include "infpar.h"
C
      PARAMETER ( NTABLE = 20, D1 = 1.0D0)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER TSTCHA*1,WORD2*6,WORD3*6
      CHARACTER*80 INLINE
#include "cbirea.h"
#include "ccom.h"
C
      DATA TABLE /'.QM3   ', '.THRDIP', '.MAXDIP', '.FIXMOM',
     *            '.OLDTG ', '.PRINT ', '.ATMVDW', '.CLONLY',
     *            '.VDWSKP', '.SKIPNC', '.MYSOLV', '.DAMPIN',
     *            '.PRFQM3', '.INTDIR', '.REDCNT', '.LGSPOL',
     *            '.MMPCM ', '.MXITMP', '.THRSMP', '.QMDAMP'/
C
      CALL QENTER('QM3INP')
C
C     Initializing QM3 keywords
C
      QM3     = .FALSE.
      THDISC  = 1.0D-9 
      MXDIIT  = 80
      FIXMOM  = .FALSE.
      OLDTG   = .FALSE.
      LONEPAR = .FALSE.
      LTWOPAR = .FALSE.
      LEPSADD = .FALSE.
      LSIGADD = .FALSE.
      LOCLAS  = .FALSE.
      IQM3PR  = IPRUSR
      VDWSKP  = .FALSE.
      SKIPNC  = .FALSE.
      MYMAT   = .FALSE.
      MYITE   = .TRUE.
      EXPON   = .FALSE.
      PRFQM3  = .FALSE.
      INTDIR  = .FALSE.
      REDCNT  = .FALSE.
      MMPCM   = .FALSE.
      MXITMP  = 20
      THRSMP  = 1.0D-05 ! MM/PCM convergence threhold
      LGSPOL  = .FALSE.
      QMDAMP  = .FALSE.
C
C
C     CRUCIAL TEST : MXCENT has to be equal to MXQM3
C
      IF (MXCENT .NE. MXQM3) THEN
        CALL QUIT('MXCENT has to be equal to MXQM3')
      END IF
C
      NEWDEF = (WORD .EQ. '*QM3   ')
      ICHANG = 0
C
      IF (NEWDEF) THEN
        WORD1 = WORD
C
  100   CONTINUE
C
        READ (LUCMD, '(A7)') WORD
        CALL UPCASE(WORD)
        PROMPT = WORD(1:1)
C
        IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
          GO TO 100
        ELSE IF (PROMPT .EQ. '.') THEN
          ICHANG = ICHANG + 1
          DO 200 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
              GO TO ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     *               11,12,13,14,15,16,17,18,19,20), I
            END IF
  200     CONTINUE
          IF (WORD .EQ. '.OPTION') THEN
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            GO TO 100
          END IF
C
          WRITE (LUPRI,'(/3A/)') ' Keyword "',WORD,
     *                           '" not recognized for *QM3'
          CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
          CALL QUIT('Illegal keyword for *QM3')
C
    1     CONTINUE
            QM3 = .TRUE.
            NYQMMM=.FALSE.
          GO TO 100
C
    2     CONTINUE
            READ(LUCMD,*) THDISC
          GO TO 100
C
    3     CONTINUE
            READ(LUCMD,*) MXDIIT
          GO TO 100
C
    4     CONTINUE
            FIXMOM = .TRUE.
          GO TO 100
C
    5     CONTINUE
            OLDTG = .TRUE.
          GO TO 100
C
    6     CONTINUE
            READ(LUCMD,*) IQM3PR
          GO TO 100
C
    7     CONTINUE 
            READ(LUCMD,'(A)') INLINE
            DO 500 K=1,80
              TSTCHA = INLINE(K:K)
              IF (TSTCHA .NE. ' ') THEN
                WORD2 = INLINE(K:K+5)
              GOTO 501
              END IF
  500       CONTINUE
  501       CONTINUE
C
          IF (WORD2 .EQ. 'TWOPAR') THEN
            LTWOPAR = .TRUE.
            GOTO 100
          ELSE
            LONEPAR = .TRUE.
            IF (WORD2 .EQ. 'EPSADD') THEN
              LEPSADD = .TRUE.
            ELSE IF (WORD2 .EQ. 'EPSMLT') THEN
              LEPSADD = .FALSE.
            ELSE IF (WORD2 .EQ. 'SIGADD') THEN
              LSIGADD = .TRUE.
            ELSE IF (WORD2 .EQ. 'SIGMLT') THEN
              LSIGADD = .FALSE.
            ELSE
              CALL QUIT('Error 1 in .ATMVDW input i DALTON.INP')
            END IF
C
            WORD3 = WORD2
C
            READ(LUCMD,'(A)') INLINE
            DO 510 K=1,80
              TSTCHA = INLINE(K:K)
              IF (TSTCHA .NE. ' ') THEN
                WORD2 = INLINE(K:K+5)
              GOTO 511
              END IF
  510       CONTINUE
  511       CONTINUE
C
            IF (WORD2 .EQ. 'EPSADD') THEN
              LEPSADD = .TRUE.
            ELSE IF (WORD2 .EQ. 'EPSMLT') THEN
              LEPSADD = .FALSE.
            ELSE IF (WORD2 .EQ. 'SIGADD') THEN
              LSIGADD = .TRUE.
            ELSE IF (WORD2 .EQ. 'SIGMLT') THEN
              LSIGADD = .FALSE.
            ELSE
              CALL QUIT('Error 2 in .ATMVDW input i DALTON.INP')
            END IF
C
            IF (WORD2(1:3) .EQ. WORD3(1:3)) THEN
            CALL QUIT('Needs a rule for both Eps and sigma in .ATMVDW!')
            END IF
          END IF
C
          GO TO 100
C
    8     CONTINUE
            LOCLAS = .TRUE. 
          GO TO 100
C
    9     CONTINUE
            VDWSKP = .TRUE.
          GO TO 100
C
   10     CONTINUE
            SKIPNC = .TRUE.
          GO TO 100
C
   11     CONTINUE
            READ(LUCMD,'(A)') INLINE
            IF (INLINE(1:5) .EQ. 'MYMAT') THEN
              MYMAT = .TRUE.
              MYITE = .FALSE.
            ELSE IF (INLINE(1:5) .EQ. 'MYITE')  THEN
              MYMAT = .FALSE.
              MYITE = .TRUE.
            ELSE
              WRITE (lupri,'(/A/A)')
     &        'Unknown my-model on input. I quit. You specified:',
     &        INLINE
              CALL QUIT('Unknown my-model on input. I quit')
            END IF
          GO TO 100
C
   12     CONTINUE
            READ(LUCMD,'(A)') INLINE
            IF (INLINE .EQ. 'EXPON') THEN
              EXPON = .TRUE.
            ELSE
              WRITE (lupri,'(/A/A)')
     &        'Only exponential damping implemented. You specified:',
     &        INLINE
              CALL QUIT('Only exponential damping implemented.')
            END IF
          GO TO 100
C
   13     CONTINUE
            PRFQM3 = .TRUE.
          GO TO 100
C
   14     CONTINUE
            INTDIR = .TRUE.
          GO TO 100

   15     CONTINUE
            REDCNT = .TRUE.
          GO TO 100

   16     CONTINUE
            LGSPOL = .TRUE.
          GO TO 100

   17     CONTINUE
            MMPCM  = .TRUE.
          GO TO 100

   18     CONTINUE
            READ(LUCMD,*) MXITMP
          GO TO 100

   19     CONTINUE
            READ(LUCMD,*) THRSMP
          GO TO 100

   20   CONTINUE
            QMDAMP = .TRUE.
            READ(LUCMD,*) ADAMP
          GO TO 100
C
        ELSE IF (PROMPT .EQ. '*') THEN
          GO TO 300
        ELSE
          WRITE (LUPRI,'(/3A/)') ' Prompt "',WORD,
     *                             '" not recognized for *QM3'
          CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
          CALL QUIT('Illegal prompt for *QM3')
        END IF
      END IF
C
  300 CONTINUE
C
      IF (ICHANG .GT. 0) THEN
        CALL HEADER('Changes of defaults for *QM3   :',0)
C
        WRITE(LUPRI,'(1X,A1,30A,A1)') '+',('-',J=1,18),'+'
        WRITE(LUPRI,'(1X,A2,A6,A3,A7,A2)')
     *         '| ',' WORD:',' | ','CHANGE:',' |'
        WRITE(LUPRI,'(1X,A1,30A,A1)') '+',('-',J=1,18),'+'
C
        IF (QM3) WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ','   QM3',' | ',QM3,' |'
C
        IF (THDISC .NE. 1.0D-9) WRITE (LUPRI,'(1X,A2,A6,A3,1P,D8.1,A2)')
     *         '| ','THDISC',' | ',THDISC,' |'
C
        IF (MXDIIT .NE. 80) WRITE (LUPRI,'(1X,A2,A6,A3,I7,A2)')
     *         '| ','MXDIIT',' | ',MXDIIT,' |'
C
        IF (FIXMOM) WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ','FIXMOM',' | ',FIXMOM,' |'
C
        IF (OLDTG) WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ',' OLDTG',' | ',OLDTG,' |'
C
          WRITE (LUPRI,'(A3,A6,A3,I7,A2)')
     *           ' | ','PRINT',' | ',IQM3PR,' |'
C
        IF (LONEPAR) THEN
          WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *           '| ','ATMVDW',' | ',LONEPAR,' |'
          IF (LEPSADD) THEN
            WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *           '| ','EPSADD',' | ',LEPSADD,' |'
          ELSE IF (.NOT. LEPSADD) THEN
            WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *           '| ','EPSMLT',' | ',(.NOT.LEPSADD),' |'
          END IF
          IF (LSIGADD) THEN
            WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *           '| ','SIGADD',' | ',LSIGADD,' |'
          ELSE IF (.NOT. LSIGADD) THEN
            WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *           '| ','SIGMLT',' | ',(.NOT.LSIGADD),' |'
          END IF
        END IF
        WRITE(LUPRI,'(A2,30A1)') ' +',('-',J=1,18),'+'
C
        WRITE(LUPRI,*) ' Settings for determination of induced dipoles:'
        IF (MYITE) THEN
          WRITE(LUPRI,*) ' Iterative Method is used'
        ELSE IF (MYMAT) THEN
          WRITE(LUPRI,*) ' Matrix inversion is used'
        END IF
        IF (EXPON) THEN
          WRITE(LUPRI,*) ' Exponential Damping is active'
        END IF
        IF (MMPCM) THEN
          WRITE(LUPRI,*) ' MM/PCM interface active'
          if (nodtot.ge.1) write(lupri,*) 'WARNING: qm3pcm is
     & not parallelized, so it will be much slower than expected!'
        END IF

        IF (NODTOT .GE. 1 .AND. .NOT. INTDIR) THEN
          WRITE(LUPRI,*) ' WARNING: QM3 .INTDIR activated since',
     *                   ' this is a parallel calculation'
          INTDIR = .TRUE.
        ENDIF

        WRITE(LUPRI,'(A2,30A1)') ' +',('-',J=1,18),'+'

        IF (QMDAMP .AND. (.NOT. INTDIR)) THEN
          CALL QUIT('QM damping only implemented for .INTDIR option')
        ENDIF
C
        IF (MMPCM .AND. .NOT. INTDIR) THEN
          WRITE(LUPRI,*) ' WARNING: QM3 .INTDIR activated since',
     *                   ' this is a QM/MM/PCM calculation'
          INTDIR = .TRUE.
        ENDIF

      END IF
C
      CALL QEXIT('QM3INP')
      RETURN
      END
C------------------------------------------------------------------------
      SUBROUTINE QM3_SETUP()
C
C     April 2009, hjaaj
C
C     Moved setup of MM info to this new subroutine.
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "qm3.h"
C
C     ***************************************************************
C     Calculating the classical interaction terms from the given
C     geometry and interaction terms:
C     ***************************************************************
C
      CALL MMFLDS()
      CALL ELNUC(ENUQM3)
C
      IF ( OLDTG ) THEN
        DO I= 1, ISYTP
          DO J = NSYSBG(I), NSYSED(I)
            DO K = 1, NCTOT
              IF (ISUBSY(K) .EQ. J) THEN
                CHAOLD(K) = CHARGE(K)
                CHARGE(K) = 0.0D0
              END IF
            END DO
          END DO
        END DO
      END IF
      RETURN
      END
C
C------------------------------------------------------------------------
C  /* Deck rdpot */
      SUBROUTINE QM3_RDPOT()
C------------------------------------------------------------------------
C
C --- General Input for MM potentials---
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "mxcent.h"
#include "abainf.h"
#include "inftap.h"
#include "qm3.h"
C
      PARAMETER (NTABLE = 24,NTAB1 = 16)
      PARAMETER (DHALF = 0.5D0, D0 = 0.0D0)
C
      CHARACTER WORD*7, WORD1*7, PROMPT*1,TSTCHA*1,
     &          TABLE(NTABLE)*7, TABLE1(NTAB1)*7
      CHARACTER*80 INLINE
      LOGICAL   LOALIS, LOALTE
      LOGICAL   LOTEST(0:MXTYPE)
C
      REAL*8 QM3SIG(0:MXTYPE,MXQM3),QM3EPS(0:MXTYPE,MXQM3)
C
      DATA TABLE  /'.TYPE  ', '.MODEL ', '.ALPISO', '.ALPTSR',
     &             '.CHARGS', '.FILINP', '.SHAWFC', '.SIGEPS',
     &             '.XXXXXX', '.XXXXXX', '.XXXXXX', '.XXXXXX',
     &             '.XXXXXX', '.XXXXXX', '.XXXXXX', '.XXXXXX',
     &             '.XXXXXX', '.XXXXXX', '.XXXXXX', '.XXXXXX',
     &             '.XXXXXX', '.XXXXXX', '.XXXXXX', '.XXXXXX'/
C
      CALL QENTER('QM3_RDPOT')
C
C     Initialization
C    
      ISYTP = 0
C
      DO 909 I = 0, MXTYPE
        ALTXX(I) = D0
        ALTXY(I) = D0
        ALTXZ(I) = D0
        ALTYY(I) = D0
        ALTYZ(I) = D0
        ALTZZ(I) = D0
        ICHRGS(I) = 0
        ISIGEPS(I)= 0
        NUALIS(I) = 0
        DISMOD(I) = .FALSE.
        RDFILE(I) = .FALSE.
        SHAWFC(I) = .FALSE.
        LOTEST(I)= .FALSE.
C
        DO 908 J = 1, MXQM3
          QM3CHG(I,J) = D0
          QM3SIG(I,J) = D0
          QM3EPS(I,J) = D0
          ALPIMM(I,J) = D0
  908   CONTINUE 
C
  909 CONTINUE   
C
      CALL TITLER('Output from MM potential input processing','*',115)
C
C     **** Find potential input *****
C
      LUPOT = -1
      CALL GPOPEN(LUPOT,'POTENTIAL.INP','OLD',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
C
      REWIND (LUPOT,IOSTAT=IOS)
C
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
C
      READ (LUPOT,'(A7)',ERR=910) WORD
C
      IF (WORD .EQ. '**SYSTP') THEN
        GOTO 930 
      ELSE
        CALL QUIT('**SYSTP has to be the first keyword' //
     &            ' in POTENTIAL.INP')
      END IF
C
  910 CONTINUE
C
      NWARN = NWARN + 1
      WRITE (LUPRI,*) 'WARNING QM3_RDPOT: error reading potential input'
C
  930 CONTINUE
C
      WORD1 = WORD
      READ (LUPOT, '(A7)') WORD
C
      IF (WORD .EQ. '.NUMMTP') THEN
        READ (LUPOT,'(I5)',ERR=940) ISYTP
      ELSE
        CALL QUIT('Illegal keyword in POTENTIAL.INP. '//
     &            '(**SYSTP has to be followed by .NUMMTP)') 
      END IF
C
      GOTO 100
C
  940 CONTINUE
C
      WRITE (LUPRI,*) 'WARNING QM3_RDPOT: error reading '//
     &                'number of system types'
      CALL QUIT('QM3_RDPOT: Wrong input for number of system types.')
C
C     ***** Process input for COMMON  /QM3SYS/  *****
C
      IF ( ISYTP .GT. MXTYPE ) THEN
         WRITE (LUPRI,'(/A,I4/A,I4/A)')
     &   ' You specified .NUMMTP under **SYSTP as',ISYTP,
     &   ' Max number of system types is MXTYPE =',MXTYPE,
     &   ' Increase MXTYPE in "include/qm3.h" and recompile'
         CALL QUIT('Too many system types, see output')
      END IF
C
  100 DO 555 J=0,ISYTP
C
        LOALIS = .FALSE.
        LOALTE = .FALSE.
C
  156   CONTINUE 
C
        READ (LUPOT, '(A7)') WORD
C
        PROMPT = WORD(1:1)
C
        IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
          GO TO 156
        ELSE IF (PROMPT .EQ. '.') THEN
C
          DO 99 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
              GO TO (101,102,103,104,105,106,107,108,109,110,111,
     *               112,113,114,115,116,117,118,119,120,121,122,
     *               123,124), I
            END IF
  99      CONTINUE
C
          IF (WORD .EQ. '.OPTION') THEN
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            GO TO 100
          END IF
C
          WRITE (LUPRI,'(/3A/)') ' Keyword ',WORD,
     &     ' in POTENTIAL.INP is not recognized in QM3_RDPOT.'
          CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
          CALL QUIT('Illegal keyword in POTENTIAL.INP 1.')
C
  101     CONTINUE
C
          IF (J .EQ. 0) THEN
            READ (LUPOT,'(A)') INLINE
            DO 344 K=1,80
              TSTCHA = INLINE(K:K) 
              IF (TSTCHA .NE. ' ') THEN
C
                READ(INLINE(K:K),*) KNUL
C
                IF (KNUL .NE. 0) THEN
                  CALL QUIT('The QM syst. has to be the first system')
                ELSE
                  NSYSBG(J) = 0
                  NSYSED(J) = 0
                  LSTART =  NSYSED(J) + 1
                  GOTO 156
                END IF 
C
              END IF
  344       CONTINUE
          ELSE
            READ (LUPOT,'(A)') INLINE
            DO 343 K= 1,80
              TSTCHA = INLINE(K:K)
              IF (TSTCHA .EQ. '-') THEN
                READ(INLINE(1:K-1),*) NSYSBG(J)
                READ(INLINE(K+1:80),*) NSYSED(J)
C
                IF (.NOT. NSYSBG(J) .LE. NSYSED(J) ) THEN
                  CALL QUIT('The MM types have been given' //
     &                      ' in wrong order')
                END IF
C
                IF ( NSYSBG(J) .NE. LSTART ) THEN
                  CALL QUIT('The MM types have to be given' //
     &                      ' in continuous order')
                ELSE
                  LSTART = NSYSED(J) + 1
                END IF
                GOTO 156
              END IF
  343       CONTINUE
            CALL QUIT('The MM types have to be seperated by a "-"')
          END IF
          GO TO 156
C
  102     CONTINUE
          READ (LUPOT,'(A)') INLINE
          DO 346 K = 1, 80
            TSTCHA = INLINE(K:K)
            IF (INLINE(K:K) .NE. ' ') THEN
              MDLWRD(J) = INLINE(K:K+6)
              GOTO 156
            END IF
  346     CONTINUE
C
  103     CONTINUE
          LOALIS = .TRUE.
          READ(LUPOT,*) NUALIS(J)
          IF (NUALIS(J) .GT. 1) THEN
            DISMOD(J) = .TRUE.
          END IF
          DO 319 MN = 1 ,NUALIS(J)
            READ (LUPOT,*) ALPIMM(J,MN)
  319     CONTINUE
          GO TO 156
C
  104     CONTINUE
          LOALTE = .TRUE.
          READ (LUPOT,*) ALTXX(J),ALTXY(J),ALTXZ(J),
     &                   ALTYY(J),ALTYZ(J),ALTZZ(J)
          GO TO 156
C
  105     CONTINUE
          READ (LUPOT,*) ICHRGS(J)
          DO 418 L = 1, ICHRGS(J)
            READ(LUPOT,*) QM3CHG(J,L)
  418     CONTINUE
          GO TO 156
C
  106     CONTINUE
          RDFILE(J) = .TRUE.
          GO TO 156
C
  107     CONTINUE
          SHAWFC(J) = .TRUE.
          GO TO 156
C
  108     CONTINUE
            IF (.NOT. LONEPAR) THEN
              CALL QUIT('Atomic vdW-par. not specified in *QM3')
            END IF
            LOTEST(J) = .TRUE.
            READ(LUPOT,*) ISIGEPS(J)
            DO 515 L = 1,ISIGEPS(J)
              READ(LUPOT,*) QM3SIG(J,L), QM3EPS(J,L)
  515       CONTINUE
          GO TO 156
C
  109     CONTINUE
          GO TO 156
C
  110     CONTINUE
          GO TO 156
C
  111     CONTINUE
          GO TO 156
C
  112     CONTINUE
          GO TO 156
C
  113     CONTINUE
          GO TO 156
C
  114     CONTINUE
          GO TO 156
C
  115     CONTINUE
          GO TO 156
C
  116     CONTINUE
          GO TO 156
C
  117     CONTINUE
          GO TO 156
C
  118     CONTINUE
          GO TO 156
C
  119     CONTINUE
          GO TO 156
C
  120     CONTINUE
          GO TO 156
C
  121     CONTINUE
          GO TO 156
C
  122     CONTINUE
          GO TO 156
C
  123     CONTINUE
          GO TO 156
C
  124     CONTINUE
          GO TO 156
C
        ELSE IF (WORD .EQ. '*******') THEN
          GO TO 444 
        ELSE
          WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal'
          CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
          CALL QUIT('Illegal prompt in POTENTIAL.INP 2.')
        END IF
C
C       Testing input according to chosen model
C
  444   IF ( (MDLWRD(J) .EQ. 'SPC    ') .AND. 
     &     ( (LOALIS) .OR. (LOALTE) ) ) THEN
          CALL QUIT('Wrong input for SPC model')
        END IF
C
        IF ( (MDLWRD(J) .EQ. 'SPC    ') .AND. (RDFILE(J)) ) THEN
          CALL QUIT('.FILINP can only be used with SPC_E models')
        END IF
C
        IF (NUALIS(J) .GT. 1) THEN
          IF ((MDLWRD(J) .EQ. 'SPC_E02') .OR.
     &        (MDLWRD(J) .EQ. 'SPC_EC2') .OR.
     &        (MDLWRD(J) .EQ. 'SPC_EC4') .OR.
     &        (MDLWRD(J) .EQ. 'SPC    ')) THEN
            CALL QUIT('No distributed alpha allowed for this model.')
          END IF
        END IF
C
        IF ( (MDLWRD(J) .EQ. 'SPC_E01') .AND.
     &    (LOALTE) .AND. (.NOT. (LOALIS)) ) THEN    
          CALL QUIT('Wrong input for SPC_E01 model 1')
        ELSE IF ( (MDLWRD(J) .EQ. 'SPC_E01') .AND.
     &    (.NOT. (LOALIS)) ) THEN
          CALL QUIT('Wrong input for SPC_E01 model 2')
        ELSE IF ( (MDLWRD(J) .EQ. 'SPC_E01') .AND.
     &    (LOALTE) .AND. (LOALIS) ) THEN
          CALL QUIT('Wrong input for SPC_E01 model 3')
        END IF
C
        IF ( (MDLWRD(J) .EQ. 'SPC_EC1') .AND.
     &    (LOALTE) .AND. (.NOT. (LOALIS)) ) THEN
          CALL QUIT('Wrong input for SPC_EC1 model 1')
        ELSE IF ( (MDLWRD(J) .EQ. 'SPC_EC1') .AND.
     &    (.NOT. (LOALIS)) ) THEN
          CALL QUIT('Wrong input for SPC_EC1 model 2')
        ELSE IF ( (MDLWRD(J) .EQ. 'SPC_EC1') .AND.
     &    (LOALTE) .AND. (LOALIS) ) THEN
          CALL QUIT('Wrong input for SPC_EC1 model 3')
        END IF
C
        IF ( (MDLWRD(J) .EQ. 'SPC_E02') .AND.
     &    (LOALIS) .AND. (.NOT. (LOALTE)) ) THEN
          CALL QUIT('Wrong input for SPC_E02 model 1')
        ELSE IF ( (MDLWRD(J) .EQ. 'SPC_E02') .AND.
     &    (.NOT. (LOALTE)) ) THEN
          CALL QUIT('Wrong input for SPC_E02 model 2')
        ELSE IF ( (MDLWRD(J) .EQ. 'SPC_E02') .AND.
     &    (LOALTE) .AND. (LOALIS) ) THEN
          CALL QUIT('Wrong input for SPC_E02 model 3')
        END IF
C
        IF ( (MDLWRD(J) .EQ. 'SPC_EC2') .AND.
     &    (LOALIS) .AND. (.NOT. (LOALTE)) ) THEN
          CALL QUIT('Wrong input for SPC_EC2 model 1')
        ELSE IF ( (MDLWRD(J) .EQ. 'SPC_EC2') .AND.
     &    (.NOT. (LOALTE)) ) THEN
          CALL QUIT('Wrong input for SPC_EC2 model 2')
        ELSE IF ( (MDLWRD(J) .EQ. 'SPC_EC2') .AND.
     &    (LOALTE) .AND. (LOALIS) ) THEN
          CALL QUIT('Wrong input for SPC_EC2 model 3')
        END IF
C
        IF ( (MDLWRD(J) .EQ. 'SPC_EC3') .AND.
     &    (LOALTE) .AND. (.NOT. (LOALIS)) ) THEN
          CALL QUIT('Wrong input for SPC_EC3 model 1')
        ELSE IF ( (MDLWRD(J) .EQ. 'SPC_EC3') .AND.
     &    (.NOT. (LOALIS)) ) THEN
          CALL QUIT('Wrong input for SPC_EC3 model 2')
        ELSE IF ( (MDLWRD(J) .EQ. 'SPC_EC3') .AND.
     &    (LOALTE) .AND. (LOALIS) ) THEN
          CALL QUIT('Wrong input for SPC_EC3 model 3')
        END IF
C
        IF ( (MDLWRD(J) .EQ. 'SPC_EC4') .AND.
     &    (LOALIS) .AND. (.NOT. (LOALTE)) ) THEN
          CALL QUIT('Wrong input for SPC_EC4 model 1')
        ELSE IF ( (MDLWRD(J) .EQ. 'SPC_EC4') .AND.
     &    (.NOT. (LOALTE)) ) THEN
          CALL QUIT('Wrong input for SPC_EC4 model 2')
        ELSE IF ( (MDLWRD(J) .EQ. 'SPC_EC4') .AND.
     &    (LOALTE) .AND. (LOALIS) ) THEN
          CALL QUIT('Wrong input for SPC_EC4 model 3')
        END IF
C
        IF (J .EQ. ISYTP) NTOTIN = NSYSED(J) 
C
C Testing for atomic sigma and epsilon vdw-parameters:
C
        IF ( (LONEPAR) .AND. (.NOT. LOTEST(J)) ) THEN
          CALL QUIT('No sigma and epsilon par. '//
     *              'given for system')
        END IF
C
  555 CONTINUE
C
C***********************************
C Interaction parameter input table:
C***********************************
C
      DATA TABLE1 /'.LJ_A  ', '.LJ_B  ', '.LJ_Aa ', '.LJ_Ba ',
     &             '.SIGEPS', 'xXXXXXX', 'xXXXXXX', 'xXXXXXX',
     &             'xXXXXXX', 'xXXXXXX', 'xXXXXXX', 'xXXXXXX',
     &             'xXXXXXX', 'xXXXXXX', 'xXXXXXX', 'xXXXXXX'/
C
C Initialization of Two-body interaction parameters:
C
      LUVDWSE = -1
      DO JBG = 0,MXTYPE
        DO IBG = 0,MXTYPE
          QM3LJA(IBG,JBG) = D0
          QM3LJB(IBG,JBG) = D0
        END DO
      END DO
C
      READ (LUPOT,'(A7)') WORD
C
      IF (WORD .EQ. '**TWOIA') THEN
        IF (LONEPAR) THEN
          CALL QUIT('No twobody interaction terms needed '//
     *    'when atomic van der Waals parameters given!')
        END IF
        GOTO 911
      ELSE IF ((WORD .EQ. '**END O') .AND.
     &         (LONEPAR)) THEN
        GOTO 666
      ELSE
        CALL QUIT('**TWOIA not found in POTENTIAL.INP')
      END IF
C
  911 CONTINUE
C
      WORD1 = WORD
      NWARN = NWARN + 1
C
      ISYTP1 = (ISYTP+1)*ISYTP/2 + ISYTP
C
  700 READ (LUPOT, '(A7)') WORD
C
      PROMPT = WORD(1:1)
C
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
        GO TO 700
      ELSE IF (PROMPT .EQ. '.') THEN
        DO 98 I = 1, NTAB1
          IF (TABLE1(I) .EQ. WORD) THEN
            GO TO (250,251,252,253,254,255,256,257,
     *             258,259,260,261,262,263,264,265), I
          END IF
  98    CONTINUE
C
        IF (WORD .EQ. '.OPTION') THEN
          CALL PRTAB(NTAB1,TABLE1,WORD1//' input keywords',LUPRI)
          GO TO 700
        END IF
C
        WRITE (LUPRI,'(/,3A,/)')
     *         ' Keyword ',WORD,' not recognized in QM3_RDPOT(2).'
        CALL PRTAB(NTAB1,TABLE1,WORD1//' input keywords',LUPRI)
        CALL QUIT('Illegal keyword in POTENTIAL.INP 3.')
C
  250   CONTINUE
        IF (LTWOPAR) THEN
          CALL QUIT('Needs atomic parameters in **TWOIA')
        END IF
C
        READ(LUPOT,*) M
C
        IF (M .NE. ISYTP1) THEN
          WRITE(LUPRI,'(A,I5)')
     *    ' Number of LJ_A parameters needed: ',ISYTP1
           CALL QUIT('Error in number of Lennard-Jones parameters A')
        END IF
C
        DO 22 I = 0,ISYTP
          DO 33 J = I,ISYTP
            IF (J .EQ. 0) GOTO 33
            READ(LUPOT,*) QM3LJA(I,J)
  33      CONTINUE
  22    CONTINUE 
        GO TO 700
C
  251   CONTINUE
        IF (LTWOPAR) THEN
          CALL QUIT('Needs atomic parameters in **TWOIA')
        END IF
C
        READ(LUPOT,*) M
C
        IF (M .NE. ISYTP1) THEN
          WRITE(LUPRI,'(A,I5)')
     *    ' Number of LJ_A parameters needed: ',ISYTP1
          CALL QUIT('Error in number of Lennard-Jones parameters B')
        END IF
C
        DO 44 I = 0,ISYTP 
          DO 55 J = I,ISYTP
            IF (J .EQ. 0) GOTO 55
            READ(LUPOT,*) QM3LJB(I,J)
  55      CONTINUE      
  44    CONTINUE
        GO TO 700
C
  252   CONTINUE
        IF (LTWOPAR) THEN
          CALL QUIT('Needs atomic parameters in **TWOIA')
        END IF
C
        READ(LUPOT,*) M
        GO TO 700
C
  253   CONTINUE
        IF (LTWOPAR) THEN
          CALL QUIT('Needs atomic parameters in **TWOIA')
        END IF
C
        READ(LUPOT,*) M
        GO TO 700
C
  254   CONTINUE
          IF (LONEPAR) THEN
            CALL QUIT('Needs no input from **TWOIA')
          END IF
          READ(LUPOT,*) NSIGEPS
          CALL GPOPEN(LUVDWSE,'VDWSIEP','UNKNOWN',' ','FORMATTED',
     &                IDUMMY,.FALSE.)
          DO 1500 KL = 1, NSIGEPS
            READ(LUPOT,*) I,K,J,L,TEMP,TEMP1
            WRITE(LUVDWSE,'(4I6,2E25.15)')
     *            I,K,J,L,TEMP,TEMP1
 1500     CONTINUE
          CALL GPCLOSE(LUVDWSE,'KEEP')
        GO TO 700
C
  255   CONTINUE
        GO TO 700
C
  256   CONTINUE
        GO TO 700
C
  257   CONTINUE
        GO TO 700
C
  258   CONTINUE
        GO TO 700
C
  259   CONTINUE
        GO TO 700
C
  260   CONTINUE
        GO TO 700
C
  261   CONTINUE
        GO TO 700
C
  262   CONTINUE
        GO TO 700
C
  263   CONTINUE
        GO TO 700
C
  264   CONTINUE
        GO TO 700
C
  265   CONTINUE
        GO TO 700
C
      ELSE IF (WORD .EQ. '**END O') THEN
        GO TO 666
      ELSE
        WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal'
        CALL PRTAB(NTAB1,TABLE1,WORD1//' input keywords',LUPRI)
        CALL QUIT('Illegal prompt in POTENTIAL.INP 4.')
      END IF
C
  666 CONTINUE
C
      CALL GPCLOSE(LUPOT,'KEEP')
C
C----------------------------------------------------
C Changing the atomic sigma and epsilon parameters to
C diatomic parameters and saving these on disc:
C----------------------------------------------------
C
      IF (LONEPAR) THEN
        LCOUNT = 0
        TEMP = D0
        TEMP1= D0
        
        CALL GPOPEN(LUVDWSE,'VDWSIEP','UNKNOWN',' ','FORMATTED',
     &              IDUMMY,.FALSE.)
        DO 1600 I = 0, ISYTP
          DO 1610 K = 1, ISIGEPS(I)
            DO 1620 J = I, ISYTP
              IF ((I .EQ. 0) .AND. (J .EQ. 0)) THEN
                GOTO 1620
              ELSE
                IF (I .EQ. J) THEN
                  M = K
                ELSE
                  M = 1
                END IF
                DO 1630 L = M, ISIGEPS(J)
                  LCOUNT = LCOUNT + 1
                  IF (LSIGADD) THEN
                    TEMP = DHALF*(QM3SIG(I,K)+QM3SIG(J,L))
                  ELSE IF (.NOT. LSIGADD) THEN
                    TEMP = SQRT(QM3SIG(I,K)*QM3SIG(J,L))
                  END IF
                  IF (LEPSADD) THEN
                    TEMP1 = DHALF*(QM3EPS(I,K)+QM3EPS(J,L))
                  ELSE IF (.NOT. LEPSADD) THEN
                    TEMP1 = SQRT(QM3EPS(I,K)*QM3EPS(J,L))
                  END IF
                  WRITE(LUVDWSE,'(4I6,2E25.15)')
     *                  I,K,J,L,TEMP,TEMP1
 1630           CONTINUE
              END IF
 1620       CONTINUE
 1610     CONTINUE
 1600   CONTINUE
        NSIGEPS = LCOUNT
        CALL GPCLOSE(LUVDWSE,'KEEP')
      END IF
C
C***************************************
C  Output design for MM potential input:
C***************************************
C
      DO 818 I = 0, ISYTP
C
        IF (I .EQ. 0) THEN
          WRITE(LUPRI,'(//,1X,72A)') ('-',J=1,72) 
          WRITE(LUPRI,'(1X,A1,A14,A1,A13,A1,A9,A1,A31,A1)')
     *                '|',' QM-sys type: ','|',
     *                ' Systems:    ','|',' Model:  ','|',
     *                ' Electric properties:           ',
     *                '|'
          WRITE(LUPRI,'(1X,72A)') ('-',J=1,72)
        END IF
C
        IF (I .EQ. 1) THEN
          WRITE(LUPRI,'(//,1X,72A)') ('-',J=1,72)
          WRITE(LUPRI,'(1X,A1,A14,A1,A13,A1,A9,A1,A31,A1)')
     *                '|',' MM-sys type: ','|',
     *                ' Systems:    ','|',' Model:  ','|',
     *                ' Electric properties:           ',
     *               '|'
          WRITE(LUPRI,'(1X,72A)') ('-',J=1,72)
        END IF
C
        IF ( (I .EQ. 0) .AND.
     &    (MDLWRD(I) .EQ. 'SPC    ') ) THEN
          WRITE(LUPRI,'(1X,A1,A6,I5,A3,
     &                  A2,A1,I4,A1,I4,A1,A2,
     &                  A8,A2,A31,A1)')
     *          '|','      ',I,'   ',
     *          '| ','[',NSYSBG(I),';',NSYSED(I),']',' |',
     *          MDLWRD(I),' |',
     *          ' No classical polarization ','|'
        ELSE IF ( (I .NE. 0) .AND.
     &    (MDLWRD(I) .EQ. 'SPC    ') ) THEN
          WRITE(LUPRI,'(1X,A1,A6,I5,A3,
     &                  A2,A1,I4,A1,I4,A1,A2,
     &                  A8,A2,A31,A1)')
     *          '|','      ',I,'   ',
     *          '| ','[',NSYSBG(I),';',NSYSED(I),']',' |',
     *          MDLWRD(I),' |',
     *          ' No polarization in this model!','|'
        ELSE IF ( (I .EQ. 0) .AND.
     &    (MDLWRD(I) .EQ. 'SPC_EC1') .OR.
     &    (MDLWRD(I) .EQ. 'SPC_EC3') ) THEN
          WRITE(LUPRI,'(1X,A1,A6,I5,A3,
     &                  A2,A1,I4,A1,I4,A1,A2,
     &                  A8,A2,A31,A1)')
     *          '|','      ',I,'   ',
     *          '| ','[',NSYSBG(I),';',NSYSED(I),']',' |',
     *          MDLWRD(I),' |',
     *          ' Charges for classical calc.:  ','|'
C
          DO 616 L = 1,ICHRGS(I)
            WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,1X,
     &                    A4,I2,A2,F8.4,14X,A1)')
     *                  '|','|','|','|','  Q(',L,')=',
     *                  QM3CHG(I,L),'|'
  616     CONTINUE
C
          IF (.NOT. (DISMOD(I))) THEN
            WRITE(LUPRI,'(1X,A1,A14,A1,A13,A1,A9,
     &                    A1,A31,A1)')
     *            '|','             ','|','            ',
     *            '|','         ','|',
     *            ' Isotropic perturb. pol.:      ','|'
          ELSE
            WRITE(LUPRI,'(1X,A1,A14,A1,A13,A1,A9,
     &                    A1,A31,A1)')
     *            '|','             ','|','            ',
     *            '|','         ','|',
     *            ' Distributed pert. pol.:       ','|'
          END IF
C
          DO 780 L = 1, NUALIS(I)
            WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,1X,
     &                    A4,I2,A2,F8.4,14X,A1)')
     *                  '|','|','|','|','alp(',L,')=',
     *                  ALPIMM(I,L),'|'
  780     CONTINUE
C
        ELSE IF ( (I .NE. 0) .AND.
     &    (MDLWRD(I) .EQ. 'SPC_EC1') .OR.
     &    (MDLWRD(I) .EQ. 'SPC_EC3') ) THEN
          IF (.NOT. (DISMOD(I))) THEN
            WRITE(LUPRI,'(1X,A1,A6,I5,A3,A2,A1,I4,
     &                    A1,I4,A1,A2,A8,A2,A31,A1)')
     *            '|','      ',I,'   ',
     *            '| ','[',NSYSBG(I),';',NSYSED(I),']',' |',
     *            MDLWRD(I),' |',
     *            ' Isotropic perturb. pol.:      ','|'
          ELSE
            WRITE(LUPRI,'(1X,A1,A6,I5,A3,A2,A1,I4,
     &                    A1,I4,A1,A2,A8,A2,A31,A1)')
     *            '|','      ',I,'   ',
     *            '| ','[',NSYSBG(I),';',NSYSED(I),']',' |',
     *            MDLWRD(I),' |',
     *            ' Distributed pert. pol.:       ','|'
          END IF
C
          DO 781 L = 1, NUALIS(I)
            WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,1X,
     &                    A4,I2,A2,F8.4,14X,A1)')
     *                  '|','|','|','|','alp(',L,')=',
     *                  ALPIMM(I,L),'|'
  781     CONTINUE
C
        ELSE IF ( (I .EQ. 0) .AND.
     &    (MDLWRD(I) .EQ. 'SPC_E01') ) THEN
          WRITE(LUPRI,'(1X,A1,A6,I5,A3,
     &                  A2,A1,I4,A1,I4,A1,A2,
     &                  A8,A2,A31,A1)')
     *          '|','      ',I,'   ',
     *          '| ','[',NSYSBG(I),';',NSYSED(I),']',' |',
     *          MDLWRD(I),' |',
     *          '                               ','|'
C
          DO 617 L = 1,ICHRGS(I)
            WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,1X,
     &                    A4,I2,A2,F8.4,14X,A1)')
     *                  '|','|','|','|','  Q(',L,')=',
     *                  QM3CHG(I,L),'|'
  617     CONTINUE
C
          IF (.NOT. (DISMOD(I))) THEN
            WRITE(LUPRI,'(1X,A1,A14,A1,A13,A1,A9,
     &                    A1,A31,A1)')
     *            '|','             ','|','            ',
     *            '|','         ','|',
     *            ' Isotropic polarizability:     ','|'
          ELSE
            WRITE(LUPRI,'(1X,A1,A14,A1,A13,A1,A9,
     &                    A1,A31,A1)')
     *            '|','             ','|','            ',
     *            '|','         ','|',
     *            ' Distributed polarizability:   ','|'
          END IF
C
          DO 782 L = 1, NUALIS(I)
            WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,1X,
     &                    A4,I2,A2,F8.4,14X,A1)')
     *                  '|','|','|','|','alp(',L,')=',
     *                  ALPIMM(I,L),'|'
  782     CONTINUE
C
        ELSE IF ( (I .NE. 0) .AND.
     &    (MDLWRD(I) .EQ. 'SPC_E01') ) THEN
          IF (.NOT. (DISMOD(I))) THEN
            WRITE(LUPRI,'(1X,A1,A6,I5,A3,A2,A1,I4,A1,I4,
     &                    A1,A2,A8,A2,A31,A1)')
     *            '|','      ',I,'   ',
     *            '| ','[',NSYSBG(I),';',NSYSED(I),']',' |',
     *            MDLWRD(I),' |',
     *            ' Isotropic polarizability:     ','|'
          ELSE
            WRITE(LUPRI,'(1X,A1,A6,I5,A3,A2,A1,I4,A1,I4,
     &                    A1,A2,A8,A2,A31,A1)')
     *            '|','      ',I,'   ',
     *            '| ','[',NSYSBG(I),';',NSYSED(I),']',' |',
     *            MDLWRD(I),' |',
     *            ' Distributed polarizability:   ','|'
          END IF
C
          DO 783 L = 1, NUALIS(I)
            WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,1X,
     &                    A4,I2,A2,F8.4,14X,A1)')
     *                  '|','|','|','|','alp(',L,')=',
     *                  ALPIMM(I,L),'|'
  783     CONTINUE
C
        ELSE IF ( (I .EQ. 0) .AND. 
     &    (MDLWRD(I) .EQ. 'SPC_EC2') .OR.
     &    (MDLWRD(I) .EQ. 'SPC_EC4') ) THEN
          WRITE(LUPRI,'(1X,A1,A6,I5,A3,
     &                  A2,A1,I4,A1,I4,A1,A2,
     &                  A8,A2,A31,A1)')
     *          '|','      ',I,'   ',
     *          '| ','[',NSYSBG(I),';',NSYSED(I),']',' |',
     *          MDLWRD(I),' |',
     *          '                               ','|'
C
          DO 618 L = 1,ICHRGS(I)
            WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,1X,
     &                    A4,I2,A2,F8.4,14X,A1)')
     *                  '|','|','|','|','  Q(',L,')=',
     *                  QM3CHG(I,L),'|'
  618     CONTINUE
C
          WRITE(LUPRI,'(1X,A1,A14,A1,A13,A1,A9,
     &                  A1,A31,A1)')
     *          '|','             ','|','            ',
     *          '|','         ','|',
     *          ' Tensor pertur. polarizability:','|'
          WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,
     &                  1X,F8.4,2X,F8.4,2X,F8.4,2X,A1)')
     *          '|','|','|','|',ALTXX(I),ALTXY(I),ALTXZ(I),'|'
          WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,
     &                  1X,F8.4,2X,F8.4,2X,F8.4,2X,A1)')
     *          '|','|','|','|',ALTXY(I),ALTYY(I),ALTYZ(I),'|'
          WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,
     &                  1X,F8.4,2X,F8.4,2X,F8.4,2X,A1)')
     *          '|','|','|','|',ALTXZ(I),ALTYZ(I),ALTZZ(I),'|'
C
        ELSE IF ( (I .NE. 0) .AND. 
     &    (MDLWRD(I) .EQ. 'SPC_EC2') .OR.
     &    (MDLWRD(I) .EQ. 'SPC_EC4') ) THEN
          WRITE(LUPRI,'(1X,A1,A6,I5,A3,A2,A1,I4,A1,I4,A1,A2,
     &                  A8,A2,A31,A1)')
     *                '|','      ',I,'   ',
     *                '| ','[',NSYSBG(I),';',NSYSED(I),']',' |',
     *                MDLWRD(I),' |',
     *                ' Tensor pertur. polarizability:','|'
          WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,
     &                  1X,F8.4,2X,F8.4,2X,F8.4,2X,A1)')
     *          '|','|','|','|',ALTXX(I),ALTXY(I),ALTXZ(I),'|'
          WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,
     &                  1X,F8.4,2X,F8.4,2X,F8.4,2X,A1)')
     *          '|','|','|','|',ALTXY(I),ALTYY(I),ALTYZ(I),'|'
          WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,
     &                  1X,F8.4,2X,F8.4,2X,F8.4,2X,A1)')
     *          '|','|','|','|',ALTXZ(I),ALTYZ(I),ALTZZ(I),'|'
C
        ELSE IF ( (I .EQ. 0) .AND. 
     &    (MDLWRD(I) .EQ. 'SPC_E02') ) THEN
          WRITE(LUPRI,'(1X,A1,A6,I5,A3,
     &                  A2,A1,I4,A1,I4,A1,A2,
     &                  A8,A2,A31,A1)')
     *          '|','      ',I,'   ',
     *          '| ','[',NSYSBG(I),';',NSYSED(I),']',' |',
     *          MDLWRD(I),' |',
     *          '                               ','|'
C
          DO 619 L = 1,ICHRGS(I)
            WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,1X,
     &                    A4,I2,A2,F8.4,14X,A1)')
     *                  '|','|','|','|','  Q(',L,')=',
     *                  QM3CHG(I,L),'|'
  619     CONTINUE
C
          WRITE(LUPRI,'(1X,A1,A14,A1,A13,A1,A9,
     &                  A1,A31,A1)')
     *          '|','             ','|','            ',
     *          '|','         ','|',
     *          ' Tensor polarizability:        ','|'
          WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,
     &                  1X,F8.4,2X,F8.4,2X,F8.4,2X,A1)')
     *          '|','|','|','|',ALTXX(I),ALTXY(I),ALTXZ(I),'|'
          WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,
     &                  1X,F8.4,2X,F8.4,2X,F8.4,2X,A1)')
     *          '|','|','|','|',ALTXY(I),ALTYY(I),ALTYZ(I),'|'
          WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,
     &                  1X,F8.4,2X,F8.4,2X,F8.4,2X,A1)')
     *          '|','|','|','|',ALTXZ(I),ALTYZ(I),ALTZZ(I),'|'
C
        ELSE IF ( (I .NE. 0) .AND. 
     &    (MDLWRD(I) .EQ. 'SPC_E02') ) THEN
          WRITE(LUPRI,'(1X,A1,A6,I5,A3,A2,A1,I4,A1,I4,A1,A2,
     &                  A8,A2,A31,A1)')
     *                '|','      ',I,'   ',
     *                '| ','[',NSYSBG(I),';',NSYSED(I),']',' |',
     *                MDLWRD(I),' |',
     *                ' Tensor polarizability:        ','|'
          WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,
     &                  1X,F8.4,2X,F8.4,2X,F8.4,2X,A1)')
     *          '|','|','|','|',ALTXX(I),ALTXY(I),ALTXZ(I),'|'
          WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,
     &                  1X,F8.4,2X,F8.4,2X,F8.4,2X,A1)')
     *          '|','|','|','|',ALTXY(I),ALTYY(I),ALTYZ(I),'|'
          WRITE(LUPRI,'(1X,A1,14X,A1,13X,A1,9X,A1,
     &                  1X,F8.4,2X,F8.4,2X,F8.4,2X,A1)')
     *          '|','|','|','|',ALTXZ(I),ALTYZ(I),ALTZZ(I),'|'
C
        ELSE
          CALL QUIT('No such MM representations implemented')
        END IF
C
        WRITE(LUPRI,'(1X,72A)') ('-',J=1,72)
C
  818 CONTINUE
C
      WRITE(LUPRI,'(/)')
C
C     Testing for any shadow wavefunctions 
C
      LG = 0
      DO 87 I=1,ISYTP
       IF (SHAWFC(I)) LG = LG +1
  87  CONTINUE
C
      IF (LG .GE. 1) THEN
        LOSHAW = .TRUE.
        WRITE(LUPRI,*)'-----------------------------------------'
        WRITE(LUPRI,*)'**** Shadow wavefunction calculation ****'
        WRITE(LUPRI,*)'-----------------------------------------'
      END IF
      CALL QEXIT('QM3_RDPOT')
      RETURN
      END
C
C----------------------------------------------------------------
C  /* Deck qm3mas */
C----------------------------------------------------------------
      SUBROUTINE QM3MAS(NONT,NONTYP)
#include "implicit.h"
C
C JK+AO nov. 2002
#include "priunit.h"
#include "mxcent.h"
#include "inftap.h"
#include "qm3.h"
C
      PARAMETER  (D0 = 0.0D0)
      INTEGER    NUALTS, KUTEST
      CHARACTER  TMPNAM*4, ID3*1
      LOGICAL    ANGQM3
      DIMENSION  NONT(*)
C
#include "nuclei.h"
#include "chrnos.h"
#include "codata.h"
C
      TOTMAS = D0
      LTEMP1 = 0
C
C-----------------------------------------------------
C     First we find the total number of molecules:
C     In this implementation it is essential that the
C     MM molecule-numbers enters the programme as
C     (n + 1) not missing any numbers in between!
C-----------------------------------------------------
C
      NTOTQM3 = 0
C
      DO 11 I = 1,NCTOT
        IF (ISUBSY(I) .GT. NTOTQM3) THEN
          NTOTQM3 = ISUBSY(I)
        END IF
  11  CONTINUE
C
C--------------------------------------------------
C     Testing if type numbers are correctly given:
C--------------------------------------------------
C
      IF (NTOTIN .NE. NTOTQM3) THEN
        CALL QUIT('Error in total no. of MM molecules')
      END IF
C
C-------------------------------------------------------
C     Testing if number of coordinates will be exceeded:
C-------------------------------------------------------
C
      KTEMP = 0
C
      DO 111 I = 0, ISYTP
        KTEMP = KTEMP + NUALIS(I) * ( NSYSED(I) - NSYSBG(I) + 1)
  111 CONTINUE
C
      IF ((KTEMP + NCTOT) .GT. MXCENT) THEN
        WRITE(LUPRI,'(//1X,A21,A11,I5)')'Number of coordinates',
     &        ' required: ',(NCTOT + KTEMP)
        WRITE(LUPRI,'(1X,A32,I5)')
     &        'Maximum number of atoms allowed:',MXCENT
        CALL QUIT('Maximum number of coordinates exceeded')
      END IF
C
C---------------------------------------------------------------
C     Checking if any coordinate input is to be read from file.
C---------------------------------------------------------------
C
      KG = 0
      DO 43 I=0,ISYTP
        IF (RDFILE(I)) KG = KG + 1
   43 CONTINUE
C
      LUOSCR = -1
      IF (KG .GT. 0) THEN
         CALL GPOPEN(LUOSCR,'QM3CORD','OLD',' ','FORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND(LUOSCR)
      END IF
C
C-------------------------------------------------------
C     Calculating the center of mass of each MM molecule
C     which does not have a distributed alpha.
C-------------------------------------------------------
C
      KURT = 0
      KURT1 = 0
C
      DO 44 L = 0, ISYTP
        IF (MDLWRD(L)(1:5) .EQ. 'SPC_E') THEN
          IF ( (.NOT. DISMOD(L)) .AND. 
     &         (.NOT. RDFILE(L)) ) THEN
            DO 45 M = NSYSBG(L), NSYSED(L)
C
              TOTALM = D0
              COMX = D0
              COMY = D0
              COMZ = D0
              LTEMP = 0
C
              DO 55 J=1,NCTOT
                IF (ISUBSY(J) .EQ. M) THEN
                  LTEMP = LTEMP + 1
C
                  CALL MASFIND(TMPMAS,J)
C
                  COMX   = COMX   + TMPMAS * CORD(1,J)
                  COMY   = COMY   + TMPMAS * CORD(2,J)
                  COMZ   = COMZ   + TMPMAS * CORD(3,J)
                  TOTALM = TOTALM + TMPMAS
                END IF
  55          CONTINUE
C
              KURT = KURT + 1
              KURT1 = KURT1 + 1
              TMINV = 1 / TOTALM
C
              CORD(1,NCTOT+KURT) = TMINV * COMX
              CORD(2,NCTOT+KURT) = TMINV * COMY
              CORD(3,NCTOT+KURT) = TMINV * COMZ
              CHARGE(NCTOT+KURT) = 0.0000
              IZATOM(NCTOT+KURT) = -2
              ISUBSY(NCTOT+KURT) = M
              ISUBSI(NCTOT+KURT) = LTEMP + 1
              LTEMP1 = LTEMP
C
              IF ((KURT - 1) .LT. 10) THEN
                NAMN(NCTOT+KURT) = 
     *               'a  '//CHRNOS(MOD((KURT1-1),10))
              ELSE IF (((KURT-1) .GT. 9) .AND. 
     &                 ((KURT-1) .LT. 100)) THEN
                NAMN(NCTOT+KURT) = 
     *               'a '//CHRNOS((KURT1-1)/10)
     &                   //CHRNOS(MOD((KURT1-1),10))
              ELSE IF (((KURT-1) .GT. 99) .AND.
     &                 ((KURT-1) .LT. 1000)) THEN
                NAMN(NCTOT+KURT) = 
     *               'a'//CHRNOS((KURT1-1)/100)
     &                  //CHRNOS(MOD((KURT1-1),100)/10)
     &                  //CHRNOS(MOD((MOD((KURT1-1),100)),10))
              ELSE
                CALL QUIT('KURT is too big (.ge.1000)')
              END IF
C
              NAMEX(3*(NCTOT+KURT)-2) = NAMN(NCTOT+KURT)//' x'
              NAMEX(3*(NCTOT+KURT)-1) = NAMN(NCTOT+KURT)//' y'
              NAMEX(3*(NCTOT+KURT))   = NAMN(NCTOT+KURT)//' z'
  45        CONTINUE
C
            NSISY(L) = LTEMP1
C
          ELSE IF (DISMOD(L) .AND. (.NOT. RDFILE(L)) ) THEN
            DO 46 M = NSYSBG(L), NSYSED(L)
C
              LTEMP = 0
C
              DO 47 J = 1, NCTOT
                IF (ISUBSY(J) .EQ. M ) THEN
C
                  LTEMP = LTEMP + 1 
                  KURT = KURT + 1
                  KURT1 = KURT1 + 1
C
                  CORD(1,NCTOT+KURT) = CORD(1,J)
                  CORD(2,NCTOT+KURT) = CORD(2,J)
                  CORD(3,NCTOT+KURT) = CORD(3,J)
                  CHARGE(NCTOT+KURT) = 0.0000
                  IZATOM(NCTOT+KURT) = -2
                  ISUBSY(NCTOT+KURT) = M
                  ISUBSI(NCTOT+KURT) = NUALIS(L) + LTEMP
C
                  IF ((KURT - 1) .LT. 10) THEN
                    NAMN(NCTOT+KURT) = 
     &                'a  '//CHRNOS(MOD((KURT1-1),10))
                  ELSE IF (((KURT-1) .GT. 9) .AND.
     &                     ((KURT-1) .LT. 100)) THEN
                    NAMN(NCTOT+KURT) = 
     &                'a '//CHRNOS((KURT1-1)/10)
     &                    //CHRNOS(MOD((KURT1-1),10))
                  ELSE IF (((KURT-1) .GT. 99) .AND.
     &                     ((KURT-1) .LT. 1000)) THEN
                    NAMN(NCTOT+KURT) = 
     &                'a'//CHRNOS((KURT1-1)/100)
     &                   //CHRNOS(MOD((KURT1-1),100)/10)
     &                   //CHRNOS(MOD((MOD((KURT1-1),100)),10))
                  ELSE
                    CALL QUIT('KURT is too big (.ge.1000)')
                  END IF
C
                  NAMEX(3*(NCTOT+KURT)-2) = NAMN(NCTOT+KURT)//' x'
                  NAMEX(3*(NCTOT+KURT)-1) = NAMN(NCTOT+KURT)//' y'
                  NAMEX(3*(NCTOT+KURT))   = NAMN(NCTOT+KURT)//' z'
                END IF
  47          CONTINUE
C
              NSISY(L) = LTEMP
C
  46        CONTINUE
C
          ELSE IF ( RDFILE(L) ) THEN
            READ (LUOSCR,'(2I3,A1)') KUTEST, NUALTS, ID3
C
C---------------------------------------------------
C           Testing if fileinput is correctly given:
C---------------------------------------------------
C
            KTEMP1 = 0  
            KTEMP1 = NUALIS(L) * ( NSYSED(L) - NSYSBG(L) + 1)
C
            IF (KUTEST .NE. L) CALL QUIT('Error in QM3CORD file')
C
            IF ( NUALTS .NE. KTEMP1 ) THEN 
              CALL QUIT('Error in QM3CORD file (no. of alpha sites).')
            END IF
C
            IF (ID3 .NE. ' ') THEN
              ANGQM3 = .TRUE.
              WRITE(LUPRI,'(/1X,A15,I4)')'Molecule type: ', L
              WRITE (LUPRI,'(/1X,A/10X,A,F11.8,A2)')
     &            'Alpha sites are entered in Angstroms'//
     &            ' and converted to atomic units.',
     &            '- Conversion factor : 1 bohr =',XTANG,' A'
            ELSE
              ANGQM3 = .FALSE.
            END IF
C
            DO 48 M = NSYSBG(L), NSYSED(L)
C
              LTEMP = 0
              LTEMP1 = 0
C
              DO 49 JL = 1, NCTOT
                IF (ISUBSY(JL) .EQ. M) THEN
                  LTEMP = LTEMP + 1
                END IF
   49         CONTINUE
C
              LTEMP1 = LTEMP
C
              DO 50 KJ = 1, NUALIS(L)
C
                LTEMP1 = LTEMP1 + 1
                KURT = KURT + 1
                KURT1 = KURT1 + 1
C
                READ (LUOSCR,*)
     *          CORD(1,NCTOT+KURT),CORD(2,NCTOT+KURT),
     *          CORD(3,NCTOT+KURT)
C      
                IF (ANGQM3) THEN
                  CORD(1,NCTOT+KURT) = CORD(1,NCTOT+KURT)/XTANG
                  CORD(2,NCTOT+KURT) = CORD(2,NCTOT+KURT)/XTANG
                  CORD(3,NCTOT+KURT) = CORD(3,NCTOT+KURT)/XTANG
                END IF
C
                CHARGE(NCTOT+KURT) = 0.0000
                IZATOM(NCTOT+KURT) = -2
                ISUBSY(NCTOT+KURT) = M
                ISUBSI(NCTOT+KURT) = LTEMP1
C
                IF ((KURT - 1) .LT. 10) THEN
                  NAMN(NCTOT+KURT) = 
     &             'a  '//CHRNOS(MOD((KURT1-1),10))
                ELSE IF (((KURT-1) .GT. 9) .AND.
     &                  ((KURT-1) .LT. 100)) THEN
                  NAMN(NCTOT+KURT) = 
     &             'a '//CHRNOS((KURT1-1)/10)
     &                 //CHRNOS(MOD((KURT1-1),10))
                ELSE IF (((KURT-1) .GT. 99) .AND.
     &                  ((KURT-1) .LT. 1000)) THEN
                  NAMN(NCTOT+KURT) = 
     &              'a'//CHRNOS((KURT1-1)/100)
     &                 //CHRNOS(MOD((KURT1-1),100)/10)
     &                 //CHRNOS(MOD((MOD((KURT1-1),100)),10))
                ELSE
                  CALL QUIT('KURT is too big (.ge.1000)')
                END IF
C
                NAMEX(3*(NCTOT+KURT)-2) = NAMN(NCTOT+KURT)//' x'
                NAMEX(3*(NCTOT+KURT)-1) = NAMN(NCTOT+KURT)//' y'
                NAMEX(3*(NCTOT+KURT))   = NAMN(NCTOT+KURT)//' z'
  50          CONTINUE
C
              NSISY(L) = LTEMP
C
  48        CONTINUE
          END IF
C
        ELSE IF (MDLWRD(L) .EQ. 'SPC    ') THEN
          DO 97 M = NSYSBG(L), NSYSED(L)
C
            LTEMP = 0
            KURT1 = KURT1 + 1
C
            DO 98 J=1,NCTOT
              IF (ISUBSY(J) .EQ. M) THEN
                LTEMP = LTEMP + 1
              END IF
  98        CONTINUE
  97      CONTINUE
C
          NSISY(L) = LTEMP
        END IF
C
  44  CONTINUE
C
      NCTOT = NCTOT + KTEMP
      NONTYP = NONTYP + 1
      NONT(NONTYP) = KTEMP
C
      IF (KG .GT. 0) THEN
        CALL GPCLOSE(LUOSCR,'KEEP')
      END IF
C
C----------------------------------------------------------
C     Testing if number of distributed polarizabilities are
C     correctly given:
C----------------------------------------------------------
C
      DO 51 L = 0, ISYTP
        IF ( DISMOD(L) .AND. (.NOT. RDFILE(L)) ) THEN 
          IF ( NSISY(L) .NE. NUALIS(L) ) THEN
            CALL QUIT('Error in no. of distributed pol. in input') 
          END IF
        END IF
   51 CONTINUE
C
C---------------------------------------------------------------------
C Testing if the number of QM nuclei is consistent with the allocation 
C in READIN in herrdn.F. This is not the most consistent way to do
C it and we have to return to this issue later on. IMPORTANT if the
C code has to be a part of the released DALTON code!!A
C---------------------------------------------------------------------
C
      IF (NSISY(0) .GT. MXCENT_QM) THEN
        WRITE(LUPRI,'(/A,I5/A,I5/A)')
     & ' Maximum number of QM atoms allowed is MXCENT_QM =',MXCENT_QM,
     & ' You have specified in input                      ',NSISY(0),
     & ' To run this calculation you must'//
     & ' increase the MXCENT_QM number in the "include/mxcent.h" file.'
        CALL QUIT(
     & 'The max. number of QM atoms allowed in mxcent.h is exceeded.')
      END IF
C
      END
C
C     ******************************************************
C  /* Deck mmflds */
      SUBROUTINE MMFLDS()
C     ******************************************************
C     This routine calculates the electric fields from the 
C     different contributers at given positions in the MM medium.
C     ******************************************************
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "qm3.h"
#include "nuclei.h"
#include "inftap.h"
C
      REAL*8     ENUC(3,MXQM3), ESITE(3,MXQM3)
      REAL*8     EQMCLF(3,MXQM3), EELEC(3,MXQM3)
      REAL*8     EIND(3,MXQM3), MYIND(3,MXQM3)
      REAL*8     EFPCM(3,MXQM3)
      PARAMETER  (D0 = 0.0D0)
C
      LOGICAL    DIPCON, LOEC3L, LOMIND 
C
C-------------------------
C     Initializing fields:
C-------------------------
C
      DO I = 1, MXQM3
        DO J = 1, 3
          ENUC(J,I)  = D0
          ESITE(J,I) = D0
          EQMCLF(J,I)= D0
          EIND(J,I)  = D0
          MYIND(J,I) = D0
          EFPCM(J,I) = D0
        END DO
      END DO
C
C     **********************************************************
C     Calculates the static electric fields 
C     due to the QM nuclei and the partial
C     charges on the MM nuclei. Also the classical field from
C     the "QM" system is calculated from the classical charges
C     given in the input file. If IPRINT .GT. 1 the electric
C     fields are printed in the output. 
C     **********************************************************
C
      NU = 0
C      
      DO 12 KI = 0, ISYTP
        IF (MDLWRD(KI)(1:5) .EQ. 'SPC_E') THEN
          DO 13 I = NSYSBG(KI), NSYSED(KI)
            DO 14 J = 1,NCTOT
              IF ( (ISUBSY(J) .EQ. I) .AND. 
     &             (ISUBSI(J) .GT. NSISY(KI)) ) THEN
C
                NU = NU + 1
C
                DO 15 K = 0, ISYTP
                  DO 21 LM = NSYSBG(K), NSYSED(K)
                    IF (LM .NE. I) THEN
                      DO 16 L = 1, NCTOT
C
                        DIST = D0
                        DISM1 = D0
                        DISM3 = D0
C
                        IF (ISUBSY(L) .EQ. LM) THEN
                          IF ( (ISUBSY(L) .EQ. 0) .AND.
     &                         (ISUBSI(L) .LE. NSISY(K)) ) THEN
C
                            DIST = (CORD(1,J)-CORD(1,L))**2 +
     &                             (CORD(2,J)-CORD(2,L))**2 + 
     &                             (CORD(3,J)-CORD(3,L))**2
                            DIST = SQRT(DIST)
                            DISM1 = 1/DIST
                            DISM3 = DISM1**3

                            IF (QMDAMP) THEN
                              DIQM = (CORD(1,J)-QMCOM(1))**2 +
     &                               (CORD(2,J)-QMCOM(2))**2 +
     &                               (CORD(3,J)-QMCOM(3))**2
                              DIQM = SQRT(DIQM)
                              DFACT = (1-exp(-ADAMP*DIQM))**3
                            ELSE
                              DFACT = 1.0D0
                            ENDIF
C
                           DO 17 M = 1, 3
                            ENUC(M,NU) = ENUC(M,NU) + 
     &                      DFACT*CHARGE(L)*(CORD(M,J)-CORD(M,L))*DISM3
                            EQMCLF(M,NU) = EQMCLF(M,NU) +
     &                      QM3CHG(ISUBSY(L),ISUBSI(L))*
     &                      (CORD(M,J)-CORD(M,L))*DISM3
  17                       CONTINUE
C
                          ELSE IF ( (ISUBSY(L) .GT. 0)
     &                      .AND. (ISUBSI(L) .LE. NSISY(K)) ) THEN 
C
                            DIST = (CORD(1,J)-CORD(1,L))**2 +
     &                             (CORD(2,J)-CORD(2,L))**2 +
     &                             (CORD(3,J)-CORD(3,L))**2
                            DIST = SQRT(DIST)
                            DISM1 = 1/DIST
                            DISM3 = DISM1**3
C                              
                            DO 18 M = 1, 3
                              ESITE(M,NU) = ESITE(M,NU) + 
     &                        CHARGE(L)*(CORD(M,J)-CORD(M,L))*DISM3
  18                        CONTINUE
C
                          END IF
                        END IF
  16                  CONTINUE
                    END IF
  21              CONTINUE
  15            CONTINUE
              END IF
  14        CONTINUE
  13      CONTINUE
        END IF
  12  CONTINUE
C
      KTEMP2 = NU
C      
C--------------------------------------
C testing for special choice of models:
C--------------------------------------
C
        LG = 0
        LH = 0
        LOEC3 = .FALSE.
        LOEC3L = .FALSE.
        LOMIND = .FALSE.
C
        DO 87 I=1,ISYTP
          IF (MDLWRD(I) .EQ. 'SPC    ') LG = LG +1
  87    CONTINUE
C
        DO 88 I=0,ISYTP
          IF (MDLWRD(I) .EQ. 'SPC_EC3') LH = LH +1
  88    CONTINUE
C
        IF (LH .GT. 0) THEN
          LOEC3 = .TRUE.
          LOEC3L = .TRUE.
          LOMIND = .FALSE.
        END IF
C
        IF (LG .EQ. ISYTP)  THEN
          WRITE(LUPRI,'(/A)') 'All MM models are SPC'
          LOSPC = .TRUE.
          IF (MDLWRD(0) .EQ. 'SPC    ') THEN
            GOTO 117
          END IF
        END IF
C
C****************************************************************
C If fixmom is true we do not want to calculate new fields
C from the induced dipole moments. In such a run we have to 
C restart from a set of known induced dipole moments. This option
C is crucial when testing the response functions as compared to
C finite field calculations.
C**************************************************************** 
C 
      IF (FIXMOM) THEN
        WRITE(LUPRI,'(1X,A38)')'Induced dipole moments are kept fixed.'
        WRITE(LUPRI,'(1X,A9,L2)')'FIXMOM = ',FIXMOM
      ELSE IF (MYITE) THEN
C
  116   CONTINUE
C
        IF (LOMIND) THEN
          LOEC3L = .FALSE.
        END IF
C
        DIPCON = .FALSE.
        XDIP1 = D0
        LM = 0
C
        DO 999 ITER = 1, MXDIIT
C
          XDIP2 = D0
          LM = LM +1
          NUM = 0
C
          DO 10 I= 0, ISYTP
            IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
              DO 20 J = NSYSBG(I), NSYSED(I)
                DO 101 K = 1, NUALIS(I)
                  NUM = NUM + 1
                  DO 40 M = 1, 3
                    IF ( .NOT. (LOMIND)) THEN
                      IF ( (MDLWRD(I) .EQ. 'SPC_E01') .OR.
     &                  (MDLWRD(I) .EQ. 'SPC_EC1') ) THEN
                        MYIND(M,NUM) = ALPIMM(I,K) *
     &                  ( EQMCLF(M,NUM) + ESITE(M,NUM) + EIND(M,NUM) )
                        XDIP2 = XDIP2 + MYIND(M,NUM)**2
                      ELSE IF (MDLWRD(I) .EQ. 'SPC_EC3') THEN
                        MYIND(M,NUM) = D0
                      END IF
                    ELSE 
                      IF (MDLWRD(I) .EQ. 'SPC_EC3') THEN
                        MYIND(M,NUM) = ALPIMM(I,K) *
     &                  ( EQMCLF(M,NUM) + ESITE(M,NUM) + EIND(M,NUM) )
                        XDIP2 = XDIP2 + MYIND(M,NUM)**2
                      END IF 
                    END IF 
  40              CONTINUE
 101            CONTINUE
C
                IF ( (MDLWRD(I) .EQ. 'SPC_E02') .OR.
     &            (MDLWRD(I) .EQ. 'SPC_EC2') .OR.
     &            (MDLWRD(I) .EQ. 'SPC_EC4') ) THEN
                  CALL QUIT('SPC_EX2 model not implemented yet')
                END IF
C
  20          CONTINUE
            END IF
  10      CONTINUE
C
          NUM = 0
C
          DO 50 I= 0, ISYTP
            IF  (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
              DO 60 J = NSYSBG(I), NSYSED(I)
                DO 70 K = 1, NCTOT
                  IF ( (ISUBSY(K) .EQ. J) .AND.
     &              (ISUBSI(K) .GT. NSISY(I)) ) THEN
C
                    NUM = NUM + 1
C
                    IF (LOMIND) THEN
                      IF (MDLWRD(I) .EQ. 'SPC_EC3') THEN
C
                        TEMPX = D0
                        TEMPY = D0
                        TEMPZ = D0
                        KURT = 0
C
                        DO 80 L = 0, ISYTP
                          IF (MDLWRD(L)(1:5) .EQ. 'SPC_E') THEN
                            DO 90 M = NSYSBG(L), NSYSED(L)
                              DO 100 N = 1, NCTOT
                                IF ( (ISUBSY(N) .EQ. M) .AND.
     &                            (ISUBSI(N) .GT. NSISY(L)) ) THEN
C
                                  KURT = KURT + 1
C
                                  IF (J .NE. M)  THEN
                                    DIST2 =(CORD(1,K)-CORD(1,N))**2
     *                                    +(CORD(2,K)-CORD(2,N))**2
     *                                    +(CORD(3,K)-CORD(3,N))**2
C
                                    DIST  = SQRT(DIST2)
                                    DIST3 = DIST2*DIST
                                    DIST4 = DIST2**2
                                    DIST5 = DIST4*DIST
C
                                    DOTPR = (CORD(1,K)-CORD(1,N))*
     &                                      MYIND(1,KURT) +
     &                                      (CORD(2,K)-CORD(2,N))*
     &                                      MYIND(2,KURT) +
     &                                      (CORD(3,K)-CORD(3,N))*
     &                                      MYIND(3,KURT)
C
                                    TEMPX = TEMPX +
     &                                      3.0*(CORD(1,K)-CORD(1,N)) *
     &                                      DOTPR/DIST5 - 
     &                                      MYIND(1,KURT)/DIST3
                                    TEMPY = TEMPY +
     &                                      3.0*(CORD(2,K)-CORD(2,N)) *
     &                                      DOTPR/DIST5 - 
     &                                      MYIND(2,KURT)/DIST3
                                    TEMPZ = TEMPZ +
     &                                      3.0*(CORD(3,K)-CORD(3,N)) *
     &                                      DOTPR/DIST5 - 
     &                                      MYIND(3,KURT)/DIST3
C
                                  END IF
                                END IF
  100                         CONTINUE
  90                        CONTINUE
                          END IF
  80                    CONTINUE
C
                        EIND(1,NUM) = TEMPX
                        EIND(2,NUM) = TEMPY
                        EIND(3,NUM) = TEMPZ
C
                      END IF
                    ELSE
                      IF (MDLWRD(I) .NE. 'SPC_EC3') THEN
C
                        TEMPX = D0
                        TEMPY = D0
                        TEMPZ = D0
                        LARS = 0
C
                        DO 81 L = 0, ISYTP
                          IF (MDLWRD(L)(1:5) .EQ. 'SPC_E') THEN
                            DO 91 M = NSYSBG(L), NSYSED(L)
                              DO 102 N = 1, NCTOT
                                IF ( (ISUBSY(N) .EQ. M) .AND.
     &                            (ISUBSI(N) .GT. NSISY(L)) ) THEN
C
                                  LARS = LARS + 1
C
                                  IF (J .NE. M) THEN
C
                                    DIST2 =(CORD(1,K)-CORD(1,N))**2 +
     *                                     (CORD(2,K)-CORD(2,N))**2 +
     *                                     (CORD(3,K)-CORD(3,N))**2
C
                                    DIST  = SQRT(DIST2)
                                    DIST3 = DIST2*DIST
                                    DIST4 = DIST2**2
                                    DIST5 = DIST4*DIST
C
                                    DOTPR = (CORD(1,K)-CORD(1,N))*
     &                                      MYIND(1,LARS) +
     &                                      (CORD(2,K)-CORD(2,N))*
     &                                      MYIND(2,LARS) +
     &                                      (CORD(3,K)-CORD(3,N))*
     &                                      MYIND(3,LARS)
C
                                    TEMPX = TEMPX +
     &                                      3.0*(CORD(1,K)-CORD(1,N)) *
     &                                      DOTPR/DIST5 - 
     &                                      MYIND(1,LARS)/DIST3
                                    TEMPY = TEMPY +
     &                                      3.0*(CORD(2,K)-CORD(2,N)) *
     &                                      DOTPR/DIST5 - 
     &                                      MYIND(2,LARS)/DIST3
                                    TEMPZ = TEMPZ +
     &                                      3.0*(CORD(3,K)-CORD(3,N)) *
     &                                      DOTPR/DIST5 - 
     &                                      MYIND(3,LARS)/DIST3
C
                                  END IF
                                END IF
  102                         CONTINUE
  91                        CONTINUE
                          END IF
  81                    CONTINUE
C
                        EIND(1,NUM) = TEMPX
                        EIND(2,NUM) = TEMPY
                        EIND(3,NUM) = TEMPZ
C
                      END IF
                    END IF
                  END IF
  70            CONTINUE
  60          CONTINUE
            END IF
  50      CONTINUE
C
          IF (ABS(XDIP1 - XDIP2) .LT. THDISC) THEN
            DIPCON = .TRUE.
            WRITE(LUPRI,'(T7,A,I4,A/T7,A,F20.15)')
     *      'QM3 induced dipole vector converged in',LM,' iterations.',
     *      'Final norm2 of QM3 induced dipole moment vector: ',XDIP2
            GOTO 1000
          ELSE
            XDIP1 = XDIP2
            IF (IQM3PR .GT. 1) WRITE(LUPRI,'(T9,A,F20.15)')
     *      'Norm2 of QM3 induced dipole moment vector:       ',XDIP2
          END IF
C
  999   CONTINUE
C
 1000   CONTINUE
C
        IF (LOEC3L) THEN
          LOMIND = .TRUE.
          WRITE(LUPRI,'(/T7,A/)')
     *    'QM3 Epol is calculated according to scheme C3'
          GOTO 116
        END IF
 
      ELSE IF (MYMAT) THEN

        IMODE = 0
        CALL MUMATINV(IMODE,KTEMP2,EQMCLF,ESITE,EFPCM,
     *                MYIND,EIND,DIPCON)

      END IF
C
      IF (.NOT. DIPCON) WRITE(LUPRI,'(/T7,A)')
     *   'QM3 induced dipole moment not converged self-consistently'
C
 117  CONTINUE
C
C     **********************************************************
C     Print section. Printing if PRINT .GT. 1.
C     **********************************************************
C
      IF (IQM3PR .GE. 2) THEN
C
C     --------------------------
C      The classical 'QM' field:
C     --------------------------
C
        WRITE(LUPRI,'(/1X,71A)') ('-',J=1,71)
        WRITE(LUPRI,'(1X,A1,A7,A1,A16,A1,A14,A1,A14,
     &  A1,A14,A1)')
     &  '|',' Field ','|',' sys, alpha site','|',
     &  '      X       ','|','      Y       ',
     &  '|','      Z       ',' |'
        WRITE(LUPRI,'(1X,71A)') ('-',J=1,71)
C
        NUM = 0
        DO 753 I = 0, ISYTP
          IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
            DO 754 J = NSYSBG(I), NSYSED(I)
              DO 755 L = 1, NUALIS(I)
                NUM = NUM + 1
                WRITE(LUPRI,'(1X,A1,A7,A1,I4,A1,I4,7X,
     &          A1,F14.10,A1,F14.10,A1,F14.10,A1)')
     &          '|',' E(QM) ','|',J,',',L,
     &          '|',EQMCLF(1,NUM),'|',EQMCLF(2,NUM),
     &          '|',EQMCLF(3,NUM),'|'
  755         CONTINUE
  754       CONTINUE
          END IF
  753   CONTINUE
C
        WRITE(LUPRI,'(1X,71A,//)') ('-',J=1,71)
C
C       ----------------
C        The site field:
C       ----------------
C
        WRITE(LUPRI,'(/1X,71A)') ('-',J=1,71)
        WRITE(LUPRI,'(A)')
     &  ' | Field | sys, alpha site|      X       |      Y       '//
     &                            '|      Z       |'
        WRITE(LUPRI,'(1X,71A)') ('-',J=1,71)
C
        NUM = 0
        DO 756 I = 0, ISYTP
          IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
            DO 757 J = NSYSBG(I), NSYSED(I)
              DO 758 L = 1, NUALIS(I)
                NUM = NUM + 1
                WRITE(LUPRI,'(1X,A1,A7,A1,I4,A1,I4,7X,
     &          A1,F14.10,A1,F14.10,A1,F14.10,A1)')
     &          '|',' Esite ','|',J,',',L,
     &          '|',ESITE(1,NUM),'|',ESITE(2,NUM),
     &          '|',ESITE(3,NUM),'|'
  758         CONTINUE 
  757       CONTINUE
          END IF
  756   CONTINUE
C
        WRITE(LUPRI,'(1X,71A,//)') ('-',J=1,71)
C
C       ----------------------
C        The nuclear QM field:
C       ----------------------
C
        WRITE(LUPRI,'(/1X,71A)') ('-',J=1,71)
        WRITE(LUPRI,'(A)')
     &  ' | Field | sys, alpha site|      X       |      Y       '//
     &                            '|      Z       |'
        WRITE(LUPRI,'(1X,71A)') ('-',J=1,71)
C
        NUM = 0 
        DO 759 I = 0, ISYTP
          IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
            DO 767 J = NSYSBG(I), NSYSED(I)
              DO 791 L = 1, NUALIS(I)
                NUM = NUM + 1
                WRITE(LUPRI,'(1X,A1,A7,A1,I4,A1,I4,7X,
     &          A1,F14.10,A1,F14.10,A1,F14.10,A1)')
     &          '|',' Enuc. ','|',J,',',L,
     &          '|',ENUC(1,NUM),'|',ENUC(2,NUM),
     &          '|',ENUC(3,NUM),'|'
  791         CONTINUE 
  767       CONTINUE
          END IF
  759   CONTINUE
C
        WRITE(LUPRI,'(1X,71A,//)') ('-',J=1,71)
C
C       ------------------------------------
C        The field from the induced dipoles:
C       ------------------------------------
C
        WRITE(LUPRI,'(/1X,71A)') ('-',J=1,71)
        WRITE(LUPRI,'(A)')
     &  ' | Field | sys, alpha site|      X       |      Y       '//
     &                            '|      Z       |'
        WRITE(LUPRI,'(1X,71A)') ('-',J=1,71)
C
        NUM = 0
        DO 761 I = 0, ISYTP
          IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
            DO 762 J = NSYSBG(I), NSYSED(I)
              DO 768 L = 1, NUALIS(I)
                NUM = NUM + 1
                WRITE(LUPRI,'(1X,A1,A7,A1,I4,A1,I4,7X,
     &          A1,F14.10,A1,F14.10,A1,F14.10,A1)')
     &          '|',' Eind  ','|',J,',',L,
     &          '|',EIND(1,NUM),'|',EIND(2,NUM),
     &          '|',EIND(3,NUM),'|'
  768         CONTINUE
  762       CONTINUE
          END IF
  761   CONTINUE
C
        WRITE(LUPRI,'(1X,71A,//)') ('-',J=1,71)
C
C       ----------------------------
C        The induced dipole moments:
C       ----------------------------
C
        WRITE(LUPRI,'(/1X,71A)') ('-',J=1,71)
        WRITE(LUPRI,'(A)')
     &  ' | Field | sys, alpha site|      X       |      Y       '//
     &                            '|      Z       |'
        WRITE(LUPRI,'(1X,71A)') ('-',J=1,71)
C
        NUM = 0
        DO 763 I = 0, ISYTP
          IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
            DO 764 J = NSYSBG(I), NSYSED(I)
              DO 769 L = 1, NUALIS(I)
                NUM = NUM + 1
                WRITE(LUPRI,'(1X,A1,A7,A1,I4,A1,I4,7X,
     &          A1,F14.10,A1,F14.10,A1,F14.10,A1)')
     &          '|',' Dipole','|',J,',',L,
     &          '|',MYIND(1,NUM),'|',MYIND(2,NUM),
     &          '|',MYIND(3,NUM),'|'
  769         CONTINUE 
  764       CONTINUE
          END IF
  763   CONTINUE
C
        WRITE(LUPRI,'(1X,71A,//)') ('-',J=1,71)
C
      END IF
C
C     **********************************************************
C     Calculating the 'classical' polarization energy:
C     **********************************************************
C
      ECLPOL = D0
      EMMPOL = D0
      NUM = 0
C
      DO 200 I = 0, ISYTP
        IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
          DO 210 J = NSYSBG(I), NSYSED(I)
            DO 220 K = 1, NUALIS(I)
              NUM = NUM + 1
                IF (I .NE. 0) THEN
                  DO 230 M = 1, 3
                    ECLPOL = ECLPOL -
     *              0.5*MYIND(M,NUM)*EQMCLF(M,NUM)
                    EMMPOL = EMMPOL -
     *              0.5*MYIND(M,NUM)*ESITE(M,NUM)
  230             CONTINUE
                END IF
  220       CONTINUE
  210     CONTINUE
        END IF
  200 CONTINUE
C
C     **********************************************************
C     Writing the calculated fields to files:
C     **********************************************************
C
C     NB NB NB NB !!!!
C
C     If FIXMOM is true non of the fields are written to
C     files. This means that both the ENUCFILE, ESITFILE and
C     the EINDFILE from a previous run must be used.
C
      LUENUC = -1
      LUSITE = -1
      LUEIND = -1
      LUPMOM = -1
      IF ( .NOT. (FIXMOM) ) THEN
        IF (MDLWRD(0) .EQ. 'SPC    ') THEN
          CALL QM3_PUT3(LUENUC,'ENUCFILE',ENUC,KTEMP2)
          CALL QM3_PUT3(LUSITE,'ESITFILE',ESITE,KTEMP2)
          CALL QM3_PUT3(LUEIND,'EINDFILE',EIND,KTEMP2)
          CALL QM3_PUT3(LUPMOM,'CLIND   ',MYIND,KTEMP2)
        ELSE
          LTEMP = MXQM3 - NUALIS(0)
          DO 911 I = 1, LTEMP
            DO 912 J = 1, 3
              ENUC(J,I)  = ENUC(J,NUALIS(0) + I)
              ESITE(J,I) = ESITE(J,NUALIS(0) + I)
              EIND(J,I)  = EIND(J,NUALIS(0) + I)
              MYIND(J,I) = MYIND(J,NUALIS(0) + I)
  912       CONTINUE
  911     CONTINUE
          DO 913 I = (LTEMP+1),MXQM3
            DO 914 J = 1,3
              ENUC(J,I)  = D0
              ESITE(J,I) = D0
              EIND(J,I)  = D0
              MYIND(J,I) = D0
  914       CONTINUE
  913     CONTINUE
          KTEMP2 = KTEMP2 - NUALIS(0)
          CALL QM3_PUT3(LUENUC,'ENUCFILE',ENUC,KTEMP2)
          CALL QM3_PUT3(LUSITE,'ESITFILE',ESITE,KTEMP2)
          CALL QM3_PUT3(LUEIND,'EINDFILE',EIND,KTEMP2)
          CALL QM3_PUT3(LUPMOM,'CLIND   ',MYIND,KTEMP2)
        END IF
      END IF
C
      END
C
C  /* Deck indmom */
      SUBROUTINE INDMOM(EELEC,LOINDM,FFJ)
C
C     Calculates the induced dipole moments by an iterative
C     procedure or by direct matrix inversion.
C
#include "implicit.h"
#include "priunit.h"
#include "qm3.h"
#include "mxcent.h"
#include "nuclei.h"
#include "inftap.h"
C
      LOGICAL    DIPCON, LOINDM 
      REAL*8     MYIND(3,MXQM3), EIND1(3,MXQM3)
      REAL*8     EIND(3,MXQM3),ESITE(3,MXQM3)
      REAL*8     ENUC(3,MXQM3),EQMCLF(3,MXQM3)
      REAL*8     EELEC(3,MXQM3),QMFIELD(3,MXQM3)
      REAL*8     FFJ(*)
      REAL*8     EFPCM(3,MXQM3)
      PARAMETER  (D0 = 0.0D0)
C
      DO 879 I = 1,MXQM3 
         DO 880 J = 1, 3
           ENUC(J,I) = D0
           ESITE(J,I) = D0
           EIND(J,I) = D0
           EIND1(J,I) = D0
           MYIND(J,I) = D0
           QMFIELD(J,I) = D0
           EFPCM(J,I) = D0
  880    CONTINUE
  879 CONTINUE
C
       IF (.NOT.(RELMOM)) PEDIP1 = 0.0D0
C
       IF (RELMOM) THEN
         WRITE(LUPRI,*)'External electric fields added to ind. mom.!!'
         IF (FFJ(1).NE.D0) WRITE(LUPRI,*)'X-Field:',FFJ(1)
         IF (FFJ(2).NE.D0) WRITE(LUPRI,*)'Y-Field:',FFJ(2)
         IF (FFJ(3).NE.D0) WRITE(LUPRI,*)'Z-Field:',FFJ(3)
       END IF
C
C     Get electric fields from files
C
      LUENUC = -1
      LUESITE = -1
      LUEIND = -1

C     In case of MM/PCM coupling calculate the PCM electric field at
C     each polarizable point

      IF (MMPCM .AND. (.NOT.LOSPC)) CALL PCMFIELD(EFPCM)

      IF (MYMAT) THEN

        CALL QM3_GET3(LUENUC,'ENUCFILE',ENUC,NCOMS)
        CALL QM3_GET3(LUESITE,'ESITFILE',ESITE,NCOMS)

        DO 881 I = 1,MXQM3
          DO 882 J = 1, 3
            QMFIELD(J,I) = ENUC(J,I) + EELEC(J,I)
  882     CONTINUE
  881   CONTINUE

C       we always just include the field from potential pcm. if not
C       present this will just be zero
 
        IMODE = 1
        CALL MUMATINV(IMODE,NCOMS,QMFIELD,ESITE,EFPCM,
     *                MYIND,EIND,DIPCON)


      ELSE IF (MYITE) THEN

        CALL QM3_GET3(LUENUC,'ENUCFILE',ENUC,NCOMS)
        CALL QM3_GET3(LUESITE,'ESITFILE',ESITE,NCOMS)
        CALL QM3_GET3(LUEIND,'EINDFILE',EIND,NCOMS)
C
        DIPCON = .FALSE.
        XDIP1 = D0
        LM = 0
        DO 999 ITER = 1, MXDIIT
          XDIP2 = D0
          LM = LM +1
C
          NUM = 0
          DO 10 I= 1, ISYTP
            IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
              DO 20 J = NSYSBG(I), NSYSED(I)
               DO 102 K = 1, NUALIS(I)
                 NUM = NUM + 1
                 DO 40 M = 1, 3
                   IF ( (MDLWRD(I) .EQ. 'SPC_E01') .OR. 
     &               (MDLWRD(I) .EQ. 'SPC_EC1') ) THEN 
                     MYIND(M,NUM) = ALPIMM(I,K) * ( EELEC(M,NUM)
     &              + ENUC(M,NUM) + ESITE(M,NUM) + EFPCM(M,NUM)
     &              + EIND(M,NUM) + FFJ(M))
                     XDIP2 = XDIP2 + MYIND(M,NUM)**2
                   ELSE IF ( (MDLWRD(I) .EQ. 'SPC_EC3') .AND.
     &               (LOINDM) )THEN
                     MYIND(M,NUM) = ALPIMM(I,K) * ( EELEC(M,NUM)
     &              + ENUC(M,NUM) + ESITE(M,NUM) + EFPCM(M,NUM) 
     &              + EIND(M,NUM) + FFJ(M))
                     XDIP2 = XDIP2 + MYIND(M,NUM)**2
                   ELSE IF ( (MDLWRD(I) .EQ. 'SPC_EC3') .AND.
     &               (.NOT. (LOINDM)) )THEN
                     MYIND(M,NUM) = D0
                   END IF
   40            CONTINUE
  102          CONTINUE
C
               IF ( (MDLWRD(I) .EQ. 'SPC_E02') .OR.
     &            (MDLWRD(I) .EQ. 'SPC_EC2') .OR.
     &            (MDLWRD(I) .EQ. 'SPC_EC4') ) THEN
                  CALL QUIT('SPC_EX2 model not implemented yet')
               END IF
  20          CONTINUE
            END IF
  10      CONTINUE
C
          NUM = 0
          DO 50 I= 1, ISYTP
            IF  (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
              DO 60 J = NSYSBG(I), NSYSED(I)
               DO 70 K = 1, NCTOT
                IF ( (ISUBSY(K) .EQ. J) .AND.
     &               (ISUBSI(K) .GT. NSISY(I)) ) THEN
                  NUM = NUM + 1
                  IF  (LOINDM) THEN  
                  IF (MDLWRD(I) .EQ. 'SPC_EC3') THEN
C
                    TEMPX = D0
                    TEMPY = D0
                    TEMPZ = D0
C
                    KURT = 0
                    DO 80 L = 1, ISYTP
                      IF (MDLWRD(L)(1:5) .EQ. 'SPC_E') THEN
                        DO 90 M = NSYSBG(L), NSYSED(L)
                          DO 100 N = 1, NCTOT
                            IF ( (ISUBSY(N) .EQ. M) .AND.
     &                           (ISUBSI(N) .GT. NSISY(L)) ) THEN
                              KURT = KURT + 1
                              IF (J .NE. M) THEN
C
                              DIST2 =((CORD(1,K)-CORD(1,N))**2 +
     *                        (CORD(2,K)-CORD(2,N))**2 +
     *                        (CORD(3,K)-CORD(3,N))**2)
C
                              DIST = SQRT(DIST2)
                              DIST3 = DIST2*DIST
                              DIST4 = DIST2**2
                              DIST5 = DIST4*DIST
C
                              DOTPR =
     &                        (CORD(1,K)-CORD(1,N))*MYIND(1,KURT) +
     &                        (CORD(2,K)-CORD(2,N))*MYIND(2,KURT) +
     &                        (CORD(3,K)-CORD(3,N))*MYIND(3,KURT)
C
                              TEMPX = TEMPX +
     &                        3.0*(CORD(1,K)-CORD(1,N)) *
     &                        DOTPR/DIST5 - MYIND(1,KURT)/DIST3
                              TEMPY = TEMPY +
     &                        3.0*(CORD(2,K)-CORD(2,N)) *
     &                        DOTPR/DIST5 - MYIND(2,KURT)/DIST3
                              TEMPZ = TEMPZ +
     &                        3.0*(CORD(3,K)-CORD(3,N)) *
     &                        DOTPR/DIST5 - MYIND(3,KURT)/DIST3
                              END IF
                            END IF
  100                     CONTINUE
  90                    CONTINUE
                      END IF
  80                CONTINUE
C
                    EIND(1,NUM) = TEMPX
                    EIND(2,NUM) = TEMPY
                    EIND(3,NUM) = TEMPZ
C
                  END IF
                  ELSE 
                    IF (MDLWRD(I) .NE. 'SPC_EC3') THEN
C
                      TEMPX = D0
                      TEMPY = D0
                      TEMPZ = D0
C
                      LARS = 0
                      DO 81 L = 1, ISYTP
                        IF (MDLWRD(L)(1:5) .EQ. 'SPC_E') THEN
                          DO 91 M = NSYSBG(L), NSYSED(L)
                            DO 101 N = 1, NCTOT
                              IF ( (ISUBSY(N) .EQ. M) .AND.
     &                             (ISUBSI(N) .GT. NSISY(L)) ) THEN
                                LARS = LARS + 1
                                IF (J .NE. M) THEN
C
                                DIST2 =((CORD(1,K)-CORD(1,N))**2 +
     *                          (CORD(2,K)-CORD(2,N))**2 +
     *                          (CORD(3,K)-CORD(3,N))**2)
C 
                                DIST = SQRT(DIST2)
                                DIST3 = DIST2*DIST
                                DIST4 = DIST2**2
                                DIST5 = DIST4*DIST
C
                                DOTPR =
     &                          (CORD(1,K)-CORD(1,N))*MYIND(1,LARS) +
     &                          (CORD(2,K)-CORD(2,N))*MYIND(2,LARS) +
     &                          (CORD(3,K)-CORD(3,N))*MYIND(3,LARS)
                                TEMPX = TEMPX +
     &                          3.0*(CORD(1,K)-CORD(1,N)) *
     &                          DOTPR/DIST5 - MYIND(1,LARS)/DIST3
                                TEMPY = TEMPY +
     &                          3.0*(CORD(2,K)-CORD(2,N)) *
     &                          DOTPR/DIST5 - MYIND(2,LARS)/DIST3
                                TEMPZ = TEMPZ +
     &                          3.0*(CORD(3,K)-CORD(3,N)) *
     &                          DOTPR/DIST5 - MYIND(3,LARS)/DIST3
                                END IF
                              END IF
  101                       CONTINUE
  91                      CONTINUE
                        END IF
  81                  CONTINUE
C
                      EIND(1,NUM) = TEMPX
                      EIND(2,NUM) = TEMPY
                      EIND(3,NUM) = TEMPZ
C
                    END IF
                  END IF
                END IF
  70           CONTINUE
  60          CONTINUE
           END IF
  50      CONTINUE
C
          IF (ABS(XDIP1 - XDIP2) .LT. THDISC) THEN
            DIPCON = .TRUE.
            WRITE(LUPRI,'(T7,A,I3,A/T7,A,F20.15)')
     *     'QM3 induced Dipole vector converged in ',LM,' iterations.',
     *     'Final norm2 of QM3 induced dipole moment vector: ',XDIP2
            GOTO 1000
          ELSE
            XDIP1 = XDIP2
            IF (IQM3PR .GT. 1) WRITE(LUPRI,'(T7,A,F20.15)')
     *       'Norm2 of QM3 induced dipole moment vector:       ',XDIP2
          END IF
C
  999 CONTINUE
 1000 CONTINUE

      END IF
C
      IF (.NOT. DIPCON) WRITE(LUPRI,'(T7,A/)')
     *   'QM3 induced dipole moment not converged self-consistently'
C
C     **********************************************************
C     Print section. Printing if PRINT .GT. 1.
C     **********************************************************
C
      IF (IQM3PR .GE. 2) THEN
C
C       ------------------------------------
C        The field from the induced dipoles:
C       ------------------------------------
C
        WRITE(LUPRI,'(/1X,71A)') ('-',J=1,71)
        WRITE(LUPRI,'(A)')
     &  ' | Field | sys, alpha site|      X       |      Y       '//
     &                            '|      Z       |'
        WRITE(LUPRI,'(1X,71A)') ('-',J=1,71)
C
        NUM = 0
        DO 761 I = 1, ISYTP
          IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
            DO 762 J = NSYSBG(I), NSYSED(I)
              DO 768 L = 1, NUALIS(I)
                NUM = NUM + 1
                WRITE(LUPRI,'(1X,A1,A7,A1,I4,A1,I4,7X,
     &          A1,F14.10,A1,F14.10,A1,F14.10,A1)')
     &          '|',' Eind  ','|',J,',',L,
     &          '|',EIND(1,NUM),'|',EIND(2,NUM),
     &          '|',EIND(3,NUM),'|'
  768         CONTINUE
  762       CONTINUE
           END IF
  761    CONTINUE
C
         WRITE(LUPRI,'(1X,71A,//)') ('-',J=1,71)
C
C       ------------------------------------
C        The total electric field:
C       ------------------------------------
C
        WRITE(LUPRI,'(/1X,71A)') ('-',J=1,71)
        WRITE(LUPRI,'(A)')
     &  ' | Field | sys, alpha site|      X       |      Y       '//
     &                            '|      Z       |'
        WRITE(LUPRI,'(1X,71A)') ('-',J=1,71)
C
        NUM = 0
        DO 361 I = 1, ISYTP
          IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
            DO 362 J = NSYSBG(I), NSYSED(I)
              DO 368 L = 1, NUALIS(I)
                NUM = NUM + 1
                WRITE(LUPRI,'(A,I4,A1,I4,7X,
     &          A1,F14.10,A1,F14.10,A1,F14.10,A1)')
     &          ' | Etot  |',J,',',L,
     &          '|',(EELEC(1,NUM)+ENUC(1,NUM)+ESITE(1,NUM)+EFPCM(1,NUM)+
     &               EIND(1,NUM)+
     &              FFJ(1)),'|',(EELEC(2,NUM)+ENUC(2,NUM)+ESITE(2,NUM)+
     &              EFPCM(2,NUM)+
     &              EIND(2,NUM)+FFJ(2)),
     &          '|',(EELEC(3,NUM)+ENUC(1,NUM)+ESITE(1,NUM)+
     &               EFPCM(3,NUM)+
     &              EIND(1,NUM)+FFJ(1)),'|'
  368         CONTINUE
  362       CONTINUE
           END IF
  361    CONTINUE
C
         WRITE(LUPRI,'(1X,71A,//)') ('-',J=1,71)
C       ----------------------------
C        The induced dipole moments:
C       ----------------------------
C
        WRITE(LUPRI,'(/1X,71A)') ('-',J=1,71)
        WRITE(LUPRI,'(A)')
     &  ' | Field | sys, alpha site|      X       |      Y       '//
     &                            '|      Z       |'
        WRITE(LUPRI,'(1X,71A)') ('-',J=1,71)
C
        NUM = 0
        DO 763 I = 1, ISYTP
          IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
            DO 764 J = NSYSBG(I), NSYSED(I)
              DO 769 L = 1, NUALIS(I)
                NUM = NUM + 1
                WRITE(LUPRI,'(1X,A1,A7,A1,I4,A1,I4,7X,
     &          A1,F14.10,A1,F14.10,A1,F14.10,A1)')
     &          '|',' Dipole','|',J,',',L,
     &          '|',MYIND(1,NUM),'|',MYIND(2,NUM),
     &          '|',MYIND(3,NUM),'|'
  769         CONTINUE 
  764       CONTINUE
           END IF
  763    CONTINUE
C
         WRITE(LUPRI,'(1X,71A)') ('-',J=1,71)
C
        END IF
C
        LUMOM = -1
        CALL QM3_PUT3(LUEIND,'EINDFILE',EIND,NCOMS)
        CALL QM3_PUT3(LUMOM,'QMIND   ',MYIND,NCOMS)
C
C       if mm/pcm we transfer the induced dipole moments to pcm        
        IF (MMPCM .AND. (.NOT.LOSPC)) CALL GIVEMY(MYIND)

        IF (SLOTH) THEN
          NUM = 0
          TEMMY1 = 0.0D0
          TEMMY2 = 0.0D0
          TEMMY3 = 0.0D0
C
          DO 263 I = 1, ISYTP
            IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
              DO 264 J = NSYSBG(I), NSYSED(I)
                DO 269 L = 1, NUALIS(I)
                  NUM = NUM + 1
                  TEMMY1 = TEMMY1 + MYIND(1,NUM)
                  TEMMY2 = TEMMY2 + MYIND(2,NUM)
                  TEMMY3 = TEMMY3 + MYIND(3,NUM)
  269           CONTINUE
  264         CONTINUE
            END IF
  263     CONTINUE
C
          WRITE(LUPRI,*)
     &    'TOTAL IND. DIPOLE MOMENT OF THE MM SYSTEM (x,y,z)'
          WRITE(LUPRI,*) TEMMY1,'au', TEMMY2,'au', TEMMY3,'au'
          WRITE(LUPRI,*)
C
          TEMPX = 0.0D0
          TEMPY = 0.0D0
          TEMPZ = 0.0D0
C
          DO 800 I = 1, ISYTP
            IF (MDLWRD(I)(1:3) .NE. 'SPC') GOTO 800
            DO 810 J = NSYSBG(I), NSYSED(I)
              DO 820 K = 1, NCTOT
                IF ( (ISUBSY(K) .EQ. J) .AND.
     *               (ISUBSI(K) .LE. NSISY(I)) ) THEN
                  TEMPX = TEMPX + CORD(1,K)*CHARGE(K)
                  TEMPY = TEMPY + CORD(2,K)*CHARGE(K)
                  TEMPZ = TEMPZ + CORD(3,K)*CHARGE(K)
                END IF
  820         CONTINUE
  810       CONTINUE
  800     CONTINUE
C
          WRITE(LUPRI,*)
     &    'TOTAL PERM. DIPOLE MOMENT OF THE MM SYSTEM (x,y,z)'
          WRITE(LUPRI,*) TEMPX,'au', TEMPY,'au', TEMPZ,'au'
          WRITE(lUPRI,*)
C
          WRITE(LUPRI,*)'TOTAL DIPOLE MOMENT OF THE MM SYSTEM (x,y,z)'
          WRITE(LUPRI,*) 
     &    TEMPX+TEMMY1,'au', TEMPY+TEMMY2,'au', TEMPZ+TEMMY3,'au'
          WRITE(lUPRI,*)
        END IF

C
        IF (RELMOM) THEN
          IF (.NOT.(SLOTH)) CALL QUIT('PUT IN .SLOTH in *CCSLV and do'//
     &    ' another run. Error from INDMOM.'//
     &    ' MM dipoles not calculated !')
          PEDIP1 = -(TEMPX*FFJ(1) + TEMPY*FFJ(2) + TEMPZ*FFJ(3))
        END IF  

      END
C    
C*******************************************************
C  /* Deck elnuc */
      SUBROUTINE ELNUC(ENUCLE)
C*******************************************************
C
C     Calculate the nuclear part of the electrostatic
C     QM/MM interaction, the classical electrostatic
C     "QM"/MM interaction, the van der Waal interaction
C     and the total classical "QM"/MM interaction.
C*******************************************************
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "qm3.h"
#include "nuclei.h"
C
      PARAMETER  (D0 = 0.0D0)
C
      EMMVDW = D0
      EMMELC = D0
      ENUCLE = D0
      ECULCL = D0
      TEMP   = D0
      TEMP1  = D0
      TEMP2  = D0
      TEMP3  = D0
      DIST   = D0
C
C******************************************************
C The classical MM/MM interaction energy is calculated:
C******************************************************
C
C The electrostatic interaction energy:
C
      DO 500 I = 1, ISYTP
        IF (MDLWRD(I)(1:3) .NE. 'SPC') GOTO 500
        DO 510 J = NSYSBG(I), NSYSED(I)
          DO 520 K = 1, NCTOT
            IF ( (ISUBSY(K) .EQ. J) .AND.
     *           (ISUBSI(K) .LE. NSISY(I)) ) THEN
              DO 530 IM = 1, ISYTP
                IF (MDLWRD(IM)(1:3) .NE. 'SPC') GOTO 530
                DO 540 JM = NSYSBG(IM), NSYSED(IM)
                  IF (JM .GT. J) THEN
                    DO 550 KM = 1, NCTOT
                      IF ( (ISUBSY(KM) .EQ. JM) .AND.
     *                     (ISUBSI(KM) .LE. NSISY(IM)) ) THEN
                        DIST = SQRT((CORD(1,K)-CORD(1,KM))**2 +
     *                        (CORD(2,K)-CORD(2,KM))**2 +
     *                        (CORD(3,K)-CORD(3,KM))**2)
                        TEMP3 = TEMP3 +
     *                          CHARGE(K)*CHARGE(KM)/DIST
                      END IF
  550               CONTINUE
                  END IF
  540           CONTINUE
  530         CONTINUE
            END IF
  520     CONTINUE
  510   CONTINUE
  500 CONTINUE
C
      EMMELC = TEMP3
C
C The van der Waals interaction energy:
C
      TEMP3 = D0
C
      IF ((LONEPAR) .OR. (LTWOPAR)) THEN
        CALL QM3ATM_VDW(TEMP3,DUMMY,1)
      ELSE
        CALL QM3VDW(TEMP3,DUMMY,1)
      END IF
C
      EMMVDW = TEMP3
C
      EMM_MM = EMMELC + EMMVDW + EMMPOL
C
C*******************************************************
C The classical "QM"/MM electrostatic energy calculated:
C*******************************************************
C
      DO 300 K=1,NCTOT
        IF ( (ISUBSY(K) .EQ. 0) .AND.
     &       (ISUBSI(K) .LE. NSISY(0)) ) THEN
          DO 310 I = 1, ISYTP
C
            IF (MDLWRD(I)(1:3) .NE. 'SPC') GOTO 310
C
            DO 320 J = NSYSBG(I), NSYSED(I)
              DO 330 L=1,NCTOT
                IF ( (ISUBSY(L) .EQ. J) .AND.
     &               (ISUBSI(L) .LE. NSISY(I)) ) THEN
                  DIST = SQRT((CORD(1,K)-CORD(1,L))**2
     *                      + (CORD(2,K)-CORD(2,L))**2
     *                      + (CORD(3,K)-CORD(3,L))**2)
                  TEMP = TEMP + CHARGE(K)*CHARGE(L)/DIST
                  TEMP1 = TEMP1 +
     *                    QM3CHG(ISUBSY(K),ISUBSI(K))*
     *                    CHARGE(L)/DIST
                END IF
  330         CONTINUE
  320       CONTINUE
  310     CONTINUE
        END IF
  300 CONTINUE
C
C The QM/MM Van der Waal energy calculated:
C
      IF ((LONEPAR) .OR. (LTWOPAR)) THEN
        CALL QM3ATM_VDW(TEMP2,EVDWSH,0)
      ELSE
        CALL QM3VDW(TEMP2,EVDWSH,0)
      END IF
C
      ENUCLE = TEMP
      ECULCL = TEMP1
      ECLVDW = TEMP2  
      ECLQM3 = ECULCL + ECLPOL + ECLVDW
C
C***********************************************************
C MM/MM energy output design:
C***********************************************************
C
      CALL TITLER('The MM/MM classical interaction energy',
     *             '*',117)
C
      WRITE(LUPRI,'(/1X,72A)') ('-',J=1,72)
      WRITE(LUPRI,'(A52,A7,F12.8,A2)')
     *' |  Eelec = Sum_n,s[ (Q_n*Q_s)/|R_n - R_s| ]        ',
     * '|      ',EMMELC,' |',
     *' |  Epol  = - 1/2*Sum_a[ Pind_a*E^site_a ]          ',
     * '|      ',EMMPOL,' |',
     *' |  Evdw  = Sum_a[ A_ma/|R_ma|^12 - B_ma/|R_ma|^6 ] ',
     * '|      ',EMMVDW,' |'
      WRITE(LUPRI,'(1X,72A)') ('-',J=1,72)
      WRITE(LUPRI,'(A52,A7,F12.8,A2)')
     *' |  E(MM/MM) = Eelec + Epol + Evdw                  ',
     * '|      ',EMM_MM,' |'
      WRITE(LUPRI,'(1X,72A//)') ('-',J=1,72)
C
C***********************************************************
C "QM"/MM energy output design:
C***********************************************************
C
      CALL TITLER('The "QM"/MM classical interaction energy',
     *             '*',115)
C
      WRITE(LUPRI,'(/1X,72A)') ('-',J=1,72)
      WRITE(LUPRI,'(A52,A7,F12.8,A2)')
     *' |  Eelec = Sum_n,s[ (Q_n*Q_s)/|R_n - R_s| ]        ',
     * '|      ',ECULCL,' |',
     *' |  Epol  = - 1/2*Sum_a[ Pind_a*E^(QMclassic)_a ]   ',
     * '|      ',ECLPOL,' |',
     *' |  Evdw  = Sum_a[ A_ma/|R_ma|^12 - B_ma/|R_ma|^6 ] ',
     * '|      ',ECLVDW,' |'
      WRITE(LUPRI,'(1X,72A)') ('-',J=1,72)
      WRITE(LUPRI,'(A52,A7,F12.8,A2)')
     *' |  E("QM"/MM) = Eelec + Epol + Evdw                ',
     * '|      ',ECLQM3,' |'
      WRITE(LUPRI,'(1X,72A//)') ('-',J=1,72)
C
C    If we have any shadow wave-functions substitute
C    ECLWDV with EVDWSH
C
      IF (LOSHAW) THEN
        ECLVDW = EVDWSH
      ENDIF
C
      IF (LOCLAS) THEN
        CALL QUIT('Exit dalton forced by keyword CLONLY for *QM3')
      END IF
C
      END
C
C*******************************************************
C  /* Deck qm3_put3*/
      SUBROUTINE QM3_PUT3(NMBU,FLNAME,PUTNU,NULOOP)
C*******************************************************
C
#include "implicit.h"
#include "mxcent.h"
#include "qm3.h"
#include "dummy.h"
#include "priunit.h"
#include "inftap.h"
C
      CHARACTER FLNAME*8
      INTEGER   NMBU,NULOOP
      DIMENSION PUTNU(3,MXQM3) 
C
      NMBU = -1
      CALL GPOPEN(NMBU,FLNAME,'UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
C
      DO 930 LM = 1,NULOOP
        WRITE(NMBU,'(I5,3E25.15)') LM,
     *  PUTNU(1,LM),PUTNU(2,LM),PUTNU(3,LM)
  930 CONTINUE
C
      CALL GPCLOSE(NMBU,'KEEP')
C
      END
C
C*******************************************************
C  /* Deck qm3_get3*/
      SUBROUTINE QM3_GET3(NMBU,FLNAME,PUTNU,NULOOP)
C*******************************************************
C
#include "implicit.h"
#include "mxcent.h"
#include "qm3.h"
#include "dummy.h"
#include "priunit.h"
#include "inftap.h"
C
      CHARACTER FLNAME*8
      INTEGER   NMBU,NULOOP
      DIMENSION PUTNU(3,MXQM3)
C
      NMBU = -1
      CALL GPOPEN(NMBU,FLNAME,'UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
C
      REWIND(NMBU)
C
      DO 920 LM = 1,NULOOP
        READ (NMBU,'(I5,3E25.15)',END=101) I,
     *        PUTNU(1,LM),PUTNU(2,LM),PUTNU(3,LM)
C
        IF (I.NE.LM) CALL QUIT( 'There is a problem in QM3_GET' )
C
        GOTO 920
  101   PUTNU(1,LM) = 1.0D0
        PUTNU(2,LM) = 1.0D0
        PUTNU(3,LM) = 1.0D0
  920 CONTINUE
C
      CALL GPCLOSE(NMBU,'KEEP')
C
      END
C
C*******************************************************
C  /* Deck qm3_obar */
      SUBROUTINE QM3_OBAR(FFJ)
C*******************************************************
C
C     Calculates the o-bar vector and the E^ns_a vector
C     entering the optimization condition of the
C     wave-function. o-bar and E^ns_a are written to 
C     files.
C
C*******************************************************
#include "implicit.h"
#include "priunit.h"
#include "qm3.h"
#include "mxcent.h"
#include "nuclei.h"
#include "inftap.h"
C
      REAL*8     EIND(3,MXQM3),ESITE(3,MXQM3)
      REAL*8     ENUC(3,MXQM3)
      REAL*8     OBAR(3,MXQM3),ENSA(3,MXQM3),EFPCM(3,MXQM3)
      REAL*8     FFJ(*)
      PARAMETER  (D0 = 0.0D0)
C
C
      DO 19 I = 1, MXQM3
         DO 20 J = 1, 3
            OBAR(J,I) = D0
            ENSA(J,I)= D0
            EFPCM(J,I) = D0
   20    CONTINUE
   19 CONTINUE
C
C---------------------------
C     get fields from files:
C---------------------------
C
C     See if external fields are to be added 
C      IF (RELMOM) THEN
C        WRITE(LUPRI,*)'External electric fields added in QM3_OBAR'
C      END IF
C
      LUENUC = -1
      LUESITE = -1
      LUEIND = -1
      LUEOBAR = -1
      LUENSA = -1
      CALL QM3_GET3(LUENUC,'ENUCFILE',ENUC,NCOMS)
      CALL QM3_GET3(LUESITE,'ESITFILE',ESITE,NCOMS)
      CALL QM3_GET3(LUEIND,'EINDFILE',EIND,NCOMS)
C
      CALL GPOPEN(LUEOBAR,'OBARFILE','UNKNOWN',' ','FORMATTED',IDUMMY,
     &           .FALSE.)
C
      CALL GPOPEN(LUENSA,'ENSAFILE','UNKNOWN',' ','FORMATTED',IDUMMY,
     &           .FALSE.)
C
      IF (MMPCM .AND. (.NOT.LOSPC)) CALL PCMFIELD(EFPCM)

C     we always just include the field from potential pcm. if not
C     present this will just be zero

      NUM = 0
C
      DO 100 I=1,ISYTP
        IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
          DO 200 J=NSYSBG(I),NSYSED(I)
            DO 300 K=1,NUALIS(I)
              NUM = NUM + 1
              DO 400 L=1,3
                OBAR(L,NUM) = 2*ENUC(L,NUM)+ESITE(L,NUM)+
     &                      EFPCM(L,NUM)+EIND(L,NUM) + FFJ(L)
                ENSA(L,NUM) = OBAR(L,NUM) + ESITE(L,NUM)+ FFJ(L) 
  400         CONTINUE
              WRITE(LUEOBAR,'(I5,3E25.15)') NUM,
     *              OBAR(1,NUM),OBAR(2,NUM),OBAR(3,NUM)
              WRITE(LUENSA,'(I5,3E25.15)') NUM,
     *              ENSA(1,NUM),ENSA(2,NUM),ENSA(3,NUM)
  300       CONTINUE
  200     CONTINUE
        END IF
  100 CONTINUE 
C
C
      CALL GPCLOSE(LUEOBAR,'KEEP')
      CALL GPCLOSE(LUENSA,'KEEP')
C
      END
C
C***************************************************
C  /* Deck qm3_otilde */
      SUBROUTINE QM3_OTILDE(OTILDE,FFJ)
C***************************************************
C
C     Calculates the o-tilde contribution to the
C     QM/MM interaction energy. 
C***************************************************
C
#include "implicit.h"
#include "priunit.h"
#include "qm3.h"
#include "mxcent.h"
#include "nuclei.h"
#include "inftap.h"
C
      REAL*8     EIND(3,MXQM3),ESITE(3,MXQM3)
      REAL*8     ENUC(3,MXQM3),RHS(3,MXQM3),EFPCM(3,MXQM3)
      REAL*8     FFJ(*)
      PARAMETER  (D0 = 0.0D0)
C
      TEMP = D0
      DOTPR = D0
      OTILDE = D0

      DO 19 I = 1, MXQM3
         DO 20 J = 1, 3
            EFPCM(J,I) = D0
   20    CONTINUE
   19 CONTINUE

C
C     See if external fields are to be added
C      IF (RELMOM) THEN
C        WRITE(LUPRI,*)'External electric fields added in QM3_OTILDE'
C      END IF  
C
C----------------------
C get field from files:
C----------------------
C
      LUENUC = -1
      LUESITE = -1
      LUEIND = -1
      CALL QM3_GET3(LUENUC,'ENUCFILE',ENUC,NCOMS)
      CALL QM3_GET3(LUESITE,'ESITFILE',ESITE,NCOMS)
      CALL QM3_GET3(LUEIND,'EINDFILE',EIND,NCOMS)
C
      IF (MMPCM .AND. (.NOT.LOSPC)) CALL PCMFIELD(EFPCM)

C     we always just include the field from potential pcm. if not
C     present this will just be zero
C
      NUM = 0
C
      DO 100 I=1,ISYTP
        IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
          DO 200 J=NSYSBG(I),NSYSED(I)
            DO 300 K = 1, NUALIS(I)
              DOTPR = D0
              NUM = NUM + 1
              DO 400 L=1,3
                RHS(L,NUM) = ENUC(L,NUM)+ESITE(L,NUM)+
     &                       EFPCM(L,NUM)+EIND(L,NUM)+FFJ(L)
                DOTPR = DOTPR + ENUC(L,NUM)*RHS(L,NUM) 
  400         CONTINUE
              OTILDE = OTILDE -0.5*ALPIMM(I,K)*DOTPR
  300       CONTINUE
  200     CONTINUE
        END IF
  100 CONTINUE
C
      END
C
C**************************************************************
C  /* Deck wrtqm3 */
      SUBROUTINE WRTQM3(AINT,LENGTH,INTTYP)
C**************************************************************
C
C     Purpose: Write MM integrals on property files.
C     Then electric field integrals are written on file
C     LUQM3E whereas the potential energy integrals are written
C     on file LUQM3P.
C     14.05.03 JK. Modified also to include overlap inegrals.
C**************************************************************
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "dummy.h"
#include "inftap.h"
#include "gnrinf.h"
#include "qm3.h"
C
      DIMENSION AINT(*)

C
      IF ( (INTTYP .EQ. 29) .AND. (QM3LO1) ) THEN
         QM3LO1 = .FALSE.
         REWIND (LUQM3E)
      ELSE IF ( (INTTYP .EQ. 35) .AND. (QM3LO2) ) THEN
         QM3LO2 = .FALSE.
         REWIND (LUQM3P)
      ELSE IF (INTTYP .EQ. 1) THEN
         REWIND (LUMMOV) 
      END IF
C
C----------------------------------
C     Write integrals to the files:
C----------------------------------
C
      IF ( INTTYP .EQ. 29) THEN
        WRITE (LUQM3E) (AINT(J), J=1,LENGTH)
      ELSE IF ( INTTYP .EQ. 35) THEN
        WRITE (LUQM3P) (AINT(J), J=1,LENGTH)
      ELSE IF ( INTTYP .EQ. 1) THEN
        WRITE (LUMMOV) (AINT(J), J=1,LENGTH)
      END IF
C
      END
C
C*************************************************************
C  /* Deck qm3_emmmm */
      SUBROUTINE QM3_EMMMM(DEMMMM,TEMPOL,FFJ)
C*************************************************************
C
#include "implicit.h"
#include "priunit.h"
#include "qm3.h"
#include "mxcent.h"
#include "nuclei.h"
#include "inftap.h"
C
      REAL*8     PIND(3,MXQM3),MYIND(3,MXQM3)
      REAL*8     ESITE(3,MXQM3),FFJ(*),EFPCM(3,MXQM3)
      PARAMETER  (D0 = 0.0D0)
C
C
C     See if external fields are to be added
C      IF (RELMOM) THEN
C       WRITE(LUPRI,*)'External electric fields added in QM3_EMMMM'
C      END IF
C
C---------------------------------------
C get induced dipole moments from files:
C---------------------------------------
C
      LUMOM = -1
      LUPMOM = -1
      LUESITE = -1
      CALL QM3_GET3(LUMOM,'QMIND   ',MYIND,NCOMS)
      CALL QM3_GET3(LUPMOM,'CLIND   ',PIND,NCOMS)
      CALL QM3_GET3(LUESITE,'ESITFILE',ESITE,NCOMS)
C
      DEMMMM = D0
      TEMPOL = D0

      DO 19 I = 1, MXQM3
         DO 20 J = 1, 3
            EFPCM(J,I) = D0
   20    CONTINUE
   19 CONTINUE

      IF (MMPCM .AND. (.NOT.LOSPC)) CALL PCMFIELD(EFPCM)

C     we always just include the field from potential pcm. if not
C     present this will just be zero

      NUM = 0
C
      DO 100 I = 1, ISYTP
        IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
          DO 110 J = NSYSBG(I), NSYSED(I)
            DO 120 L= 1, NUALIS(I)
              NUM = NUM + 1
              DO 130 K = 1, 3
                DEMMMM = DEMMMM - 0.5D0*(MYIND(K,NUM)-PIND(K,NUM))*
     *                   (ESITE(K,NUM)+EFPCM(K,NUM)+FFJ(K))
                TEMPOL = TEMPOL - 0.5D0*MYIND(K,NUM)*(ESITE(K,NUM)+
     *                   EFPCM(K,NUM)+FFJ(K))
  130         CONTINUE
  120       CONTINUE
  110     CONTINUE
        END IF
  100 CONTINUE
C
      END
C
C**************************************************************
C  /* Deck qm3vdw */
      SUBROUTINE QM3VDW(TEMP2,TEMP3,LBEGIN)
C**************************************************************
C
C     Calculate the van der Waal interaction between the QM
C     system and the MM systems when LBEGIN=0 and the Van der
C     Waal interaction between the MM systems when LBEGIN = 1.
C
C**************************************************************
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "qm3.h"
#include "nuclei.h"
C
      PARAMETER  (D0 = 0.0D0)
C
      TEMP2 = D0
      TEMP3 = D0
      TOTALM = D0
      TOTALM2 = D0
C
      IF (VDWSKP) RETURN
C
      DO 200 I = LBEGIN, ISYTP
        DO 210 J = NSYSBG(I), NSYSED(I)
          IF (DISMOD(I) .OR. (MDLWRD(I) .EQ. 'SPC    ') ) THEN
C
            QMCORX = D0
            QMCORY = D0
            QMCORZ = D0
            TOTALM = D0
C
            DO 220 K = 1, NCTOT
C
              QMMASS = D0
C
              IF ( (ISUBSY(K) .EQ. J) .AND.
     &             (ISUBSI(K) .LE. NSISY(I)) ) THEN
                CALL MASFIND(QMMASS,K)
                QMCORX = QMCORX + QMMASS * CORD(1,K)
                QMCORY = QMCORY + QMMASS * CORD(2,K)
                QMCORZ = QMCORZ + QMMASS * CORD(3,K)
                TOTALM = TOTALM + QMMASS
              END IF
  220       CONTINUE
C
            QMCORX = QMCORX/TOTALM
            QMCORY = QMCORY/TOTALM
            QMCORZ = QMCORZ/TOTALM
C
          ELSE
C
            QMCORX = D0
            QMCORY = D0
            QMCORZ = D0
C
            DO 230 K = 1, NCTOT
              IF ( (ISUBSY(K) .EQ. J) .AND.
     &             (ISUBSI(K) .GT. NSISY(I)) ) THEN
                QMCORX = CORD(1,K)
                QMCORY = CORD(2,K)
                QMCORZ = CORD(3,K)
              END IF
  230       CONTINUE
          END IF
C
          DO 240 IM = 1, ISYTP
            DO 250 JM = NSYSBG(IM), NSYSED(IM)
              IF (JM .GT. J) THEN
C
                DIST = D0
                DISM1  = D0
                DISM6  = D0
                DISM12 = D0
                TMCORX = D0
                TMCORY = D0
                TMCORZ = D0
C
                IF (DISMOD(IM) .OR. 
     &            (MDLWRD(IM) .EQ. 'SPC    ')) THEN
C
                  TOTALM2 = D0
C
                  DO 260 KM = 1, NCTOT
                    IF ( (ISUBSY(KM) .EQ. JM) .AND.
     &                   (ISUBSI(KM) .LE. NSISY(IM)) ) THEN
C
                      CALL MASFIND(TMMASS,KM)
                  
                      TMCORX = TMCORX + TMMASS * CORD(1,KM)
                      TMCORY = TMCORY + TMMASS * CORD(2,KM)
                      TMCORZ = TMCORZ + TMMASS * CORD(3,KM)
                      TOTALM2 = TOTALM2 + TMMASS
                    END IF
  260             CONTINUE
C
                  TMCORX = TMCORX/TOTALM2
                  TMCORY = TMCORY/TOTALM2
                  TMCORZ = TMCORZ/TOTALM2
C
                ELSE IF (.NOT. (DISMOD(IM)) ) THEN
                  DO 270 KM = 1, NCTOT
                    IF ( (ISUBSY(KM) .EQ. JM) .AND.
     &                   (ISUBSI(KM) .GT. NSISY(IM)) ) THEN
                      TMCORX = CORD(1,KM)
                      TMCORY = CORD(2,KM)
                      TMCORZ = CORD(3,KM)
                    END IF
  270             CONTINUE
                END IF
C
                DIST = SQRT((QMCORX-TMCORX)**2 +
     *                 (QMCORY-TMCORY)**2 +
     *                 (QMCORZ-TMCORZ)**2)
                DISTM1 = 1/DIST
                DISTM6 = DISTM1**6
                DSTM12 = DISTM6**2
                TEMP2 = TEMP2 + QM3LJA(I,IM)*DSTM12 -
     *                  QM3LJB(I,IM)*DISTM6
                IF (I .EQ. 0) THEN
                  IF (SHAWFC(IM)) THEN
                    TEMP3 = TEMP3 
                  ELSE
                    TEMP3 = TEMP3 + QM3LJA(I,IM)*DSTM12 -
     *                      QM3LJB(I,IM)*DISTM6
                  END IF
                END IF
Cftnchek        TEMJA1 = TEMJA1 + QM3LJA(I,IM)*DSTM12
Cftnchek        TEMJA2 = TEMJA2 + QM3LJB(I,IM)*DISTM6
              END IF
  250       CONTINUE
  240     CONTINUE
          IF (LBEGIN .EQ. 0) GOTO 999
  210   CONTINUE
  200 CONTINUE
  999 CONTINUE
C
      END
C
C**************************************************************
C  /* Deck qm3atm_vdw */
      SUBROUTINE QM3ATM_VDW(TEMP2,TEMP3,LBEGIN)
C**************************************************************
C
C     Calculate the van der Waal interaction between the QM
C     system and the MM systems when LBEGIN=0 and the Van der
C     Waal interaction between the MM systems when LBEGIN = 1.
C     Contrary to the QM3VDW ROUTINE this routine calculates
C     the van der Waal interaction from atomic interactions
C     according to:
C     Evdw = Sum(i<j)[ 4*Eps(ij)*[ (Sig(ij)/R(ij))**12 -
C                                  (Sig(ij)/R(ij))**6 ]]
C
C**************************************************************
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "qm3.h"
#include "nuclei.h"
C
      PARAMETER  (D0 = 0.0D0)
C
      SIG_IJ = D0
      EPS_IJ = D0
C
      LUVDWSE = -1
      CALL GPOPEN(LUVDWSE,'VDWSIEP','UNKNOWN',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND (LUVDWSE)
C
      DO 100 IVDW = 1, NSIGEPS
        READ(LUVDWSE,'(4I6,2E25.15)')
     *       I,K,J,L,SIG_IJ,EPS_IJ
        IF (LBEGIN .NE. 0) THEN
          IF (I .EQ. 0) GOTO 100
        ELSE IF (LBEGIN .EQ. 0) THEN
          IF (I .NE. 0) GOTO 100
        END IF
        DO 110 JVDW = 1, NCTOT
          IF ( (ISUBSY(JVDW) .GE. NSYSBG(I)) .AND.
     *         (ISUBSY(JVDW) .LE. NSYSED(I)) .AND.
     *         (ISUBSI(JVDW) .EQ. K) ) THEN
            DO 120 KVDW = 1, NCTOT
C
            DIST = D0
            SIG_RM1 = D0
            SIG_RM6 = D0
            SIG_RM12= D0
C
              IF ( (ISUBSY(KVDW) .GE. NSYSBG(J)) .AND.
     *             (ISUBSY(KVDW) .LE. NSYSED(J)) .AND.
     *             (ISUBSY(KVDW) .GT. ISUBSY(JVDW)) .AND.
     *             (ISUBSI(KVDW) .EQ. L) ) THEN
C
                DIST = (CORD(1,JVDW)-CORD(1,KVDW))**2 +
     *                 (CORD(2,JVDW)-CORD(2,KVDW))**2 +
     *                 (CORD(3,JVDW)-CORD(3,KVDW))**2
                DIST = SQRT(DIST)
                SIG_RM1 = SIG_IJ/DIST
                SIG_RM6 = SIG_RM1**6
                SIG_RM12= SIG_RM6**2
C
                TEMP2 = TEMP2 +
     *                  4.0D0*EPS_IJ*(SIG_RM12-SIG_RM6)
                IF (I .EQ. 0) THEN
                  IF (SHAWFC(J)) THEN
                    TEMP3 = TEMP3
                  ELSE
                    TEMP3 = TEMP3 +
     *                      4.0D0*EPS_IJ*(SIG_RM12-SIG_RM6)
                  END IF
                END IF
Cftnchek        TEMJA1 = TEMJA1 + 4.0D0*EPS_IJ*SIG_RM12
Cftnchek        TEMJA2 = TEMJA2 + 4.0D0*EPS_IJ*SIG_RM6
C
              END IF
 120        CONTINUE
          END IF
 110    CONTINUE
 100  CONTINUE
C
      IF (LBEGIN .NE. 0) THEN
        CALL GPCLOSE(LUVDWSE,'KEEP')
      ELSE IF (LBEGIN .EQ. 0) THEN
        CALL GPCLOSE(LUVDWSE,'DELETE')
      END IF
C
      END
C
C*************************************************************
C  /* Deck masfind */
      SUBROUTINE MASFIND(TMPMAS,J)
C*************************************************************
C
C     This routine finds the mass from the name of the atom
C     given in the MOLECULE.INP file.
C
C*************************************************************
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
C
      PARAMETER  (D0 = 0.0D0)
C
      IF (NAMN(J)(1:2) .EQ. 'H ') THEN
        TMPMAS = DISOTP(1,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'He') THEN
        TMPMAS = DISOTP(2,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Li') THEN
        TMPMAS = DISOTP(3,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Be') THEN
        TMPMAS = DISOTP(4,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'B ') THEN
        TMPMAS = DISOTP(5,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'C ') THEN
        TMPMAS = DISOTP(6,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'N ') THEN
        TMPMAS = DISOTP(7,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'O ') THEN
        TMPMAS = DISOTP(8,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'F ') THEN
        TMPMAS = DISOTP(9,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Ne') THEN
        TMPMAS = DISOTP(10,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Na') THEN
        TMPMAS = DISOTP(11,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Mg') THEN
        TMPMAS = DISOTP(12,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Al') THEN
        TMPMAS = DISOTP(13,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Si') THEN
        TMPMAS = DISOTP(14,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'P ') THEN
        TMPMAS = DISOTP(15,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'S ') THEN
        TMPMAS = DISOTP(16,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Cl') THEN
        TMPMAS = DISOTP(17,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Ar') THEN
        TMPMAS = DISOTP(18,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'K ') THEN
        TMPMAS = DISOTP(19,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Ca') THEN
        TMPMAS = DISOTP(20,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Sc') THEN
        TMPMAS = DISOTP(21,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Ti') THEN
        TMPMAS = DISOTP(22,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'V ') THEN
        TMPMAS = DISOTP(23,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Cr') THEN
        TMPMAS = DISOTP(24,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Mn') THEN
        TMPMAS = DISOTP(25,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Fe') THEN
        TMPMAS = DISOTP(26,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Co') THEN
        TMPMAS = DISOTP(27,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Ni') THEN
        TMPMAS = DISOTP(28,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Cu') THEN
        TMPMAS = DISOTP(29,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Zn') THEN
        TMPMAS = DISOTP(30,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Ga') THEN
        TMPMAS = DISOTP(31,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Ge') THEN
        TMPMAS = DISOTP(32,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'As') THEN
        TMPMAS = DISOTP(33,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Se') THEN
        TMPMAS = DISOTP(34,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Br') THEN
        TMPMAS = DISOTP(35,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Kr') THEN
        TMPMAS = DISOTP(36,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Rb') THEN
        TMPMAS = DISOTP(37,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Sr') THEN
        TMPMAS = DISOTP(38,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Y ') THEN
        TMPMAS = DISOTP(39,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Zr') THEN
        TMPMAS = DISOTP(40,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Nb') THEN
        TMPMAS = DISOTP(41,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Mo') THEN
        TMPMAS = DISOTP(42,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Tc') THEN
        TMPMAS = DISOTP(43,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Ru') THEN
        TMPMAS = DISOTP(44,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Rh') THEN
        TMPMAS = DISOTP(45,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Pd') THEN
        TMPMAS = DISOTP(46,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Ag') THEN
        TMPMAS = DISOTP(47,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Cd') THEN
        TMPMAS = DISOTP(48,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'In') THEN
        TMPMAS = DISOTP(49,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Sn') THEN
        TMPMAS = DISOTP(50,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Sb') THEN
        TMPMAS = DISOTP(51,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Te') THEN
        TMPMAS = DISOTP(52,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'I ') THEN
        TMPMAS = DISOTP(53,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Xe') THEN
        TMPMAS = DISOTP(54,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Cs') THEN
        TMPMAS = DISOTP(55,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Ba') THEN
        TMPMAS = DISOTP(56,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'La') THEN
        TMPMAS = DISOTP(57,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Ce') THEN
        TMPMAS = DISOTP(58,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Pr') THEN
        TMPMAS = DISOTP(59,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Nd') THEN
        TMPMAS = DISOTP(60,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Pm') THEN
        TMPMAS = DISOTP(61,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Sm') THEN
        TMPMAS = DISOTP(62,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Eu') THEN
        TMPMAS = DISOTP(63,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Gd') THEN
        TMPMAS = DISOTP(64,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Tb') THEN
        TMPMAS = DISOTP(65,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Dy') THEN
        TMPMAS = DISOTP(66,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Ho') THEN
        TMPMAS = DISOTP(67,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Er') THEN
        TMPMAS = DISOTP(68,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Tm') THEN
        TMPMAS = DISOTP(69,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Yb') THEN
        TMPMAS = DISOTP(70,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Lu') THEN
        TMPMAS = DISOTP(71,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Hf') THEN
        TMPMAS = DISOTP(72,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Ta') THEN
        TMPMAS = DISOTP(73,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'W ') THEN
        TMPMAS = DISOTP(74,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Re') THEN
        TMPMAS = DISOTP(75,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Os') THEN
        TMPMAS = DISOTP(76,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Ir') THEN
        TMPMAS = DISOTP(77,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Pt') THEN
        TMPMAS = DISOTP(78,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Au') THEN
        TMPMAS = DISOTP(79,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Hg') THEN
        TMPMAS = DISOTP(80,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Tl') THEN
        TMPMAS = DISOTP(81,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Pb') THEN
        TMPMAS = DISOTP(82,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Bi') THEN
        TMPMAS = DISOTP(83,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Po') THEN
        TMPMAS = DISOTP(84,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'At') THEN
        TMPMAS = DISOTP(85,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'Rn') THEN
        TMPMAS = DISOTP(86,1,'MASS')
      ELSE IF (NAMN(J)(1:2) .EQ. 'XX') THEN
        TMPMAS = D0
      ELSE
        WRITE(LUPRI,*)
     &    'MASFIND: Check atom no.',J,' in input file: ',NAMN(J)
        CALL QUIT('MASFIND: This atom not recognized.')
      END IF
C
      END
C
C*************************************************************
C  /* Deck qm3qcc1 */
      SUBROUTINE QM3QCC1(LUPRI,IPRINT)
C
C     Nuclear contribution to electric field gradients 
C     from the MM system. 1. contribution from the partial charges.
C
#include "implicit.h"
#include "mxcent.h"
#include "qm3.h"
      PARAMETER (D0 = 0.0D0, D3 = 3.0D0)
#include "nuclei.h"
#include "nqcc.h"
C
      IF (OLDTG)  THEN
        DO 873 I= 1, ISYTP
         DO 874 J = NSYSBG(I), NSYSED(I)
           DO 875 K = 1, NCTOT
              IF (ISUBSY(K) .EQ. J) CHARGE(K) = CHAOLD(K)
  875      CONTINUE
  874    CONTINUE
  873   CONTINUE
      END IF
C
      NATOM = 1
      DO 30 I = 1, NCTOT
        IF ( (ISUBSY(I) .EQ. 0) .AND.
     &       (ISUBSI(I) .LE. NSISY(0)) )THEN
           DO 40 I2 =1, ISYTP
           IF (MDLWRD(I2)(1:3) .NE. 'SPC') GOTO 40 
             DO 50 J = NSYSBG(I2), NSYSED(I2) 
               DO 60 L=1,NCTOT
                 IF ( (ISUBSY(L) .EQ. J) .AND.
     &                (ISUBSI(L) .LE. NSISY(I2)) ) THEN
                   IF (CHARGE(L) .EQ. 0.0D0) GOTO 60 
                   XCOOR = CORD(1,I) - CORD(1,L)
                   YCOOR = CORD(2,I) - CORD(2,L)
                   ZCOOR = CORD(3,I) - CORD(3,L)
                   R2 = XCOOR*XCOOR + YCOOR*YCOOR + ZCOOR*ZCOOR
                   R5 = R2*R2*SQRT(R2)
                   UCNNQC(1,1,NATOM) = UCNNQC(1,1,NATOM)
     &                     + CHARGE(L)*(D3*XCOOR*XCOOR - R2)/R5
                   UCNNQC(2,2,NATOM) = UCNNQC(2,2,NATOM)
     &                     + CHARGE(L)*(D3*YCOOR*YCOOR - R2)/R5
                   UCNNQC(3,3,NATOM) = UCNNQC(3,3,NATOM)
     &                     + CHARGE(L)*(D3*ZCOOR*ZCOOR - R2)/R5
                   UCNNQC(1,2,NATOM) = UCNNQC(1,2,NATOM)
     &                     + CHARGE(L)*D3*XCOOR*YCOOR/R5
                   UCNNQC(1,3,NATOM) = UCNNQC(1,3,NATOM)
     &                     + CHARGE(L)*D3*XCOOR*ZCOOR/R5
                   UCNNQC(2,3,NATOM) = UCNNQC(2,3,NATOM)
     &                     + CHARGE(L)*D3*YCOOR*ZCOOR/R5
                 END IF
 60            CONTINUE
 50          CONTINUE
 40        CONTINUE
         UCNNQC(2,1,NATOM) = UCNNQC(1,2,NATOM)
         UCNNQC(3,1,NATOM) = UCNNQC(1,3,NATOM)
         UCNNQC(3,2,NATOM) = UCNNQC(2,3,NATOM)
         NATOM = NATOM + 1
        END IF
 30   CONTINUE
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('electric field gradient after QM/MM charges',
     &               -1)
         CALL EFGPRI(UCNNQC,'                                 ')
      END IF
C
C     Finally, we set the charges to zero if OLDTG
C
      IF (OLDTG) THEN
       DO 876 I= 1, ISYTP
          DO 877 J = NSYSBG(I), NSYSED(I)
            DO 878 K = 1, NCTOT
               IF (ISUBSY(K) .EQ. J) CHARGE(K) = D0
  878       CONTINUE
  877    CONTINUE
  876   CONTINUE
      END IF
      RETURN
C
      END
C*************************************************************
C  /* Deck qm3qcc2 */
      SUBROUTINE QM3QCC2(LUPRI,IPRINT)
C
C     Nuclear contribution to electric field gradients
C     from the MM system. 2. Contribution from the induced dipoles.
C
#include "implicit.h"
#include "mxcent.h"
#include "qm3.h"
#include "nuclei.h"
#include "nqcc.h"
C
      PARAMETER (D0 = 0.0D0, D3 = 3.0D0)
      REAL*8 MYIND(3,MXQM3)
C
C
C     get induced dipole moments from files:
C
      LUMOM = -1
      CALL QM3_GET3(LUMOM,'QMIND   ',MYIND,NCOMS)
C
      IF (OLDTG)  THEN
        DO 873 I= 1, ISYTP
         DO 874 J = NSYSBG(I), NSYSED(I)
           DO 875 K = 1, NCTOT
              IF (ISUBSY(K) .EQ. J) CHARGE(K) = CHAOLD(K)
  875      CONTINUE
  874    CONTINUE
  873   CONTINUE
      END IF
C
      NATOM = 1
      DO 30 I = 1, NCTOT
        IF ( (ISUBSY(I) .EQ. 0) .AND.
     &       (ISUBSI(I) .LE. NSISY(0)) )THEN
           NUM = 0
           DO 40 I2 =1, ISYTP
           IF (MDLWRD(I2)(1:5) .NE. 'SPC_E') GOTO 40
             DO 50 J = NSYSBG(I2), NSYSED(I2)
               DO 60 L=1,NCTOT
                 IF ( (ISUBSY(L) .EQ. J) .AND.
     &                (ISUBSI(L) .GT. NSISY(I2)) ) THEN
                   TEXXX = 0.0D0
                   TEXXY = 0.0D0
                   TEXXZ = 0.0D0
                   TEYYX = 0.0D0
                   TEYYY = 0.0D0
                   TEYYZ = 0.0D0
                   TEZZX = 0.0D0
                   TEZZY = 0.0D0
                   TEZZZ = 0.0D0
                   TEXYZ = 0.0D0
                   NUM = NUM +1
                   XCOOR = CORD(1,I) - CORD(1,L)
                   YCOOR = CORD(2,I) - CORD(2,L)
                   ZCOOR = CORD(3,I) - CORD(3,L)
                   R2 = XCOOR*XCOOR + YCOOR*YCOOR + ZCOOR*ZCOOR
                   R5 = R2*R2*SQRT(R2)
                   R7 = R5*R2
                   TEXXX = -(15*XCOOR*XCOOR*XCOOR-D3*R2*D3*XCOOR)/R7
                   TEXXY = -(15*XCOOR*XCOOR*YCOOR-D3*R2*YCOOR)/R7
                   TEXXZ = -(15*XCOOR*XCOOR*ZCOOR-D3*R2*ZCOOR)/R7
                   TEYYX = -(15*YCOOR*YCOOR*XCOOR-D3*R2*XCOOR)/R7
                   TEYYY = -(15*YCOOR*YCOOR*YCOOR-D3*R2*D3*YCOOR)/R7
                   TEYYZ = -(15*YCOOR*YCOOR*ZCOOR-D3*R2*ZCOOR)/R7
                   TEZZX = -(15*ZCOOR*ZCOOR*XCOOR-D3*R2*XCOOR)/R7
                   TEZZY = -(15*ZCOOR*ZCOOR*YCOOR-D3*R2*YCOOR)/R7
                   TEZZZ = -(15*ZCOOR*ZCOOR*ZCOOR-D3*R2*D3*ZCOOR)/R7
                   TEXYZ = -(15*XCOOR*YCOOR*ZCOOR)/R7
                   UCNNQC(1,1,NATOM) = UCNNQC(1,1,NATOM)
     &                        - TEXXX*MYIND(1,NUM) - TEXXY*MYIND(2,NUM)
     &                        - TEXXZ*MYIND(3,NUM)
                   UCNNQC(2,2,NATOM) = UCNNQC(2,2,NATOM)
     &                        - TEYYX*MYIND(1,NUM) - TEYYY*MYIND(2,NUM)
     &                        - TEYYZ*MYIND(3,NUM)
                   UCNNQC(3,3,NATOM) = UCNNQC(3,3,NATOM)
     &                        - TEZZX*MYIND(1,NUM) - TEZZY*MYIND(2,NUM)
     &                        - TEZZZ*MYIND(3,NUM)
                   UCNNQC(1,2,NATOM) = UCNNQC(1,2,NATOM)
     &                        - TEXXY*MYIND(1,NUM) - TEYYX*MYIND(2,NUM)
     &                        - TEXYZ*MYIND(3,NUM)
                   UCNNQC(1,3,NATOM) = UCNNQC(1,3,NATOM)
     &                        - TEXXZ*MYIND(1,NUM) - TEXYZ*MYIND(2,NUM)
     &                        - TEZZX*MYIND(3,NUM)
                   UCNNQC(2,3,NATOM) = UCNNQC(2,3,NATOM)
     &                        - TEXYZ*MYIND(1,NUM) - TEYYZ*MYIND(2,NUM)
     &                        - TEZZY*MYIND(3,NUM)
                 END IF
 60            CONTINUE
 50          CONTINUE
 40        CONTINUE
         UCNNQC(2,1,NATOM) = UCNNQC(1,2,NATOM)
         UCNNQC(3,1,NATOM) = UCNNQC(1,3,NATOM)
         UCNNQC(3,2,NATOM) = UCNNQC(2,3,NATOM)
         NATOM = NATOM + 1
        END IF
 30   CONTINUE
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('electric field gradient after QM/MM ind. dipoles',
     &               -1)
         CALL EFGPRI(UCNNQC,'                                 ')
      END IF
C
C     Finally, we set the charges to zero if OLDTG
C
      IF (OLDTG) THEN
       DO 876 I= 1, ISYTP
          DO 877 J = NSYSBG(I), NSYSED(I)
            DO 878 K = 1, NCTOT
               IF (ISUBSY(K) .EQ. J) CHARGE(K) = D0
  878       CONTINUE
  877    CONTINUE
  876   CONTINUE
      END IF
      RETURN
C
      END
C*******************************************************
      SUBROUTINE EQMNTEST()
C*******************************************************
C
C     Calculate the QM nuclear part of the electric field
C     at specific points. If the field is requested at 
C     at specific nuclei the contribution due to that nuclei
C     is omitted.
C
C*******************************************************
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "qm3.h"
#include "nuclei.h"
C
      PARAMETER  (D0 = 0.0D0)
C
      DIST   = D0
      DIST2  = D0
      DIST3  = D0
C
      DO 100 K = 1, NCTOT
        TEMPX  = D0
        TEMPY  = D0
        TEMPZ  = D0
        DO 101 KM = 1, NCTOT
          IF (K .NE. KM) THEN 
            DIST = SQRT((CORD(1,K)-CORD(1,KM))**2 +
     *             (CORD(2,K)-CORD(2,KM))**2 +
     *             (CORD(3,K)-CORD(3,KM))**2)
            DIST2 = DIST**2
            DIST3 = DIST2*DIST
            TEMPX = TEMPX +
     *              CHARGE(KM)*( (CORD(1,K)-CORD(1,KM)) )/DIST3
            TEMPY = TEMPY +
     *              CHARGE(KM)*( (CORD(2,K)-CORD(2,KM)) )/DIST3
            TEMPZ = TEMPZ +
     *              CHARGE(KM)*( (CORD(3,K)-CORD(3,KM)) )/DIST3
          END IF 
  101   CONTINUE
            WRITE(LUPRI,*)
            WRITE(LUPRI,*)'Nuclear electric field for atom/site no.',K
            WRITE(LUPRI,*) TEMPX, TEMPY, TEMPZ
            WRITE(LUPRI,*)
  100 CONTINUE
C
      RETURN 
      END
C**********************************************************************
C  /* Deck rdqm3 */
      SUBROUTINE RDQM3(AINT,LOFRST,INTTYP)
C**************************************************************
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "dummy.h"
#include "inftap.h"
#include "inforb.h"
#include "gnrinf.h"
#include "qm3.h"
C
      DIMENSION AINT(*)
      LOGICAL LOFRST

C
      IF ( (INTTYP .EQ. 29) .AND. (LOFRST) ) THEN
         LOFRST = .FALSE.
         REWIND (LUQM3E)
      ELSE IF ( (INTTYP .EQ. 35) .AND. (LOFRST) ) THEN
         LOFRST = .FALSE.
         REWIND (LUQM3P)
      END IF
C
C----------------------------------
C     Read integrals from files:
C----------------------------------
C
      IF ( INTTYP .EQ. 29) THEN
        READ (LUQM3E) (AINT(J), J=1,NNBASX)
      ELSE IF ( INTTYP .EQ. 35) THEN
        READ (LUQM3P) (AINT(J), J=1,NNBASX)
      END IF
C
      END
C*********************************************************************
C  /* Deck mumatinv */
      SUBROUTINE MUMATINV(IMODE,KTEMP2,EQMCLF,ESITE,EFPCM,MYIND,
     &                    EIND,DIPCON)
C**************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "qm3.h"
#include "nuclei.h"
#include "inftap.h"

      INTEGER    IMODE, KTEMP2  

      REAL*8     TMAT(3*KTEMP2,3*KTEMP2), TMATI(3*KTEMP2,3*KTEMP2)
      REAL*8     EQMCLT(3*KTEMP2), TMUM(3*KTEMP2), EINDT(3*KTEMP2)
      REAL*8     EQMCLF(3,MXQM3), ESITE(3,MXQM3), MYIND(3,MXQM3) 
      REAL*8     EIND(3,MXQM3), EFPCM(3,MXQM3)

      PARAMETER  (D0 = 0.0D0)

      LOGICAL    DIPCON

C initialization

      DO 001 I=1,3*KTEMP2
        TMUM(I)  = 0.0D0
        EINDT(I) = 0.0D0
  001 CONTINUE

      KSIT1 = 1
      KSIT2 = 1

        DO 133 L = IMODE, ISYTP
          IF ( NUALIS(L) .GT. 1 ) THEN
            CALL QUIT('Matrix inversion not impl. for distributed pol')
          END IF
          IF (MDLWRD(L)(1:5) .EQ. 'SPC_E') THEN
            DO 134 M = NSYSBG(L), NSYSED(L)
              DO 135 N = 1, NCTOT
                IF ( (ISUBSY(N) .EQ. M) .AND.
     &          (ISUBSI(N) .GT. NSISY(L)) ) THEN

                 KSIT2 = 1

                  DO 136 L1 = IMODE, ISYTP
                    IF (MDLWRD(L1)(1:5) .EQ. 'SPC_E') THEN
                      DO 137 M1 = NSYSBG(L1), NSYSED(L1)
                        DO 138 N1 = 1, NCTOT
                          IF ( (ISUBSY(N1) .EQ. M1) .AND.
     &                    (ISUBSI(N1) .GT. NSISY(L1)) ) THEN

                            IF (M .EQ. M1) THEN

                              TMAT(KSIT1,KSIT2)     = -1/ALPIMM(L1,1)
                              TMAT(KSIT1,KSIT2+1)   = 0.0D0
                              TMAT(KSIT1,KSIT2+2)   = 0.0D0
                              TMAT(KSIT1+1,KSIT2)   = 0.0D0
                              TMAT(KSIT1+1,KSIT2+1) = -1/ALPIMM(L1,1)
                              TMAT(KSIT1+1,KSIT2+2) = 0.0D0
                              TMAT(KSIT1+2,KSIT2)   = 0.0D0
                              TMAT(KSIT1+2,KSIT2+1) = 0.0D0
                              TMAT(KSIT1+2,KSIT2+2) = -1/ALPIMM(L1,1)

                            ELSE

                              DIST2 =((CORD(1,N1)-CORD(1,N))**2 +
     *                               (CORD(2,N1)-CORD(2,N))**2 +
     *                               (CORD(3,N1)-CORD(3,N))**2)

                              DIST = SQRT(DIST2)
                              DIST3 = DIST2*DIST
                              DIST4 = DIST2**2
                              DIST5 = DIST4*DIST

                              IF (EXPON) THEN
                                ALPEFF1 = (ALPIMM(L,1)*ALPIMM(L1,1))
                                ALPEFF  = ALPEFF1**(1.0D0/6.0D0)
                                SCREEN  = 2.1304*DIST/ALPEFF
                                FEPQ    = 1- (1+SCREEN+0.5D0*SCREEN**2)*
     *                                    EXP(-SCREEN)
                                FTPQ    = FEPQ - (1.0D0/6.0D0)*
     *                                    SCREEN**3*EXP(-SCREEN)
                              ELSE
                                FEPQ = 1.0D0
                                FTPQ = 1.0D0
                              END IF 

                              TMAT(KSIT1,KSIT2) = 
     &                        3*FTPQ*(CORD(1,N1)-CORD(1,N))**2/DIST5 -
     &                        FEPQ*1.0D0/DIST3

                              TMAT(KSIT1,KSIT2+1) = 
     &                        3*FTPQ*(CORD(1,N1)-CORD(1,N))*
     &                        (CORD(2,N1)-CORD(2,N))/DIST5

                              TMAT(KSIT1,KSIT2+2) = 
     &                        3*FTPQ*(CORD(1,N1)-CORD(1,N))*
     &                        (CORD(3,N1)-CORD(3,N))/DIST5

                              TMAT(KSIT1+1,KSIT2) = 
     &                        3*FTPQ*(CORD(2,N1)-CORD(2,N))*
     &                        (CORD(1,N1)-CORD(1,N))/DIST5

                              TMAT(KSIT1+1,KSIT2+1) = 
     &                        3*FTPQ*(CORD(2,N1)-CORD(2,N))**2/DIST5 -
     &                        FEPQ*1.0D0/DIST3

                              TMAT(KSIT1+1,KSIT2+2) = 
     &                        3*FTPQ*(CORD(2,N1)-CORD(2,N))*
     &                        (CORD(3,N1)-CORD(3,N))/DIST5

                              TMAT(KSIT1+2,KSIT2) = 
     &                        3*FTPQ*(CORD(3,N1)-CORD(3,N))*
     &                        (CORD(1,N1)-CORD(1,N))/DIST5

                              TMAT(KSIT1+2,KSIT2+1) = 
     &                        3*FTPQ*(CORD(3,N1)-CORD(3,N))*
     &                        (CORD(2,N1)-CORD(2,N))/DIST5

                              TMAT(KSIT1+2,KSIT2+2) = 
     &                        3*FTPQ*(CORD(3,N1)-CORD(3,N))**2/DIST5 -
     &                        FEPQ*1.0D0/DIST3

                            END IF

                          KSIT2 = KSIT2 + 3

                          END IF
  138                   CONTINUE
  137                 CONTINUE
                    END IF
  136             CONTINUE

                  KSIT1 = KSIT1 + 3

                END IF
  135         CONTINUE
  134       CONTINUE
          END IF
  133   CONTINUE

      CALL DF_MATINV(TMAT,TMATI,3*KTEMP2)    

C       we always just include the field from potential pcm (EFPCM). if not
C       present this will just be zero

      JJ = 1
      DO 557 I=1,KTEMP2
        EQMCLT(JJ)   = EQMCLF(1,I) + ESITE(1,I) + EFPCM(1,I)
        EQMCLT(JJ+1) = EQMCLF(2,I) + ESITE(2,I) + EFPCM(2,I)
        EQMCLT(JJ+2) = EQMCLF(3,I) + ESITE(3,I) + EFPCM(3,I)
        JJ=JJ+3
  557 CONTINUE


      DO 558 I=1,3*KTEMP2
        DO 559 J=1,3*KTEMP2
         TMUM(I) = TMUM(I) - TMATI(I,J)*EQMCLT(J)
  559   CONTINUE
  558 CONTINUE

      I = 1
      DO 560 J=1,KTEMP2
            MYIND(1,J) = TMUM(I)
            MYIND(2,J) = TMUM(I+1)
            MYIND(3,J) = TMUM(I+2)
            I = I +3
  560 CONTINUE
C
        DO 561 J=1,3*KTEMP2
          DO 562 I=1,3*KTEMP2
           IF (J .NE. I) EINDT(J) = EINDT(J) + TMAT(J,I)*TMUM(I)
  562     CONTINUE
  561   CONTINUE

      I = 1
      DO 563 J=1,KTEMP2
            EIND(1,J) = EINDT(I)
            EIND(2,J) = EINDT(I+1)
            EIND(3,J) = EINDT(I+2)
            I = I +3
  563 CONTINUE

      DIPCON = .TRUE.
      END
C*********************************************************************
      SUBROUTINE RFQMMM()
C
C     Calculate the reaction field (RF) at the QM nuclei 
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "qm3.h"
#include "nuclei.h"
C
      PARAMETER (D0 = 0.0D0, D3 = 3.0D0)
      REAL*8     MYIND(3,MXQM3)
      REAL*8     RFQM(3,MXCENT_QM)
      CHARACTER  LAPRPC*8
C

      IF (.NOT. PRFQM3) RETURN

      DO 11 I=1,MXCENT_QM
        RFQM(1,I)=D0
        RFQM(2,I)=D0
        RFQM(3,I)=D0
  11  CONTINUE

C     if oldtg then put in charges 

      IF (OLDTG)  THEN
        DO 873 I= 1, ISYTP
         DO 874 J = NSYSBG(I), NSYSED(I)
           DO 875 K = 1, NCTOT
              IF (ISUBSY(K) .EQ. J) CHARGE(K) = CHAOLD(K)
  875      CONTINUE
  874    CONTINUE
  873   CONTINUE
      END IF
  
C     now, find a QM nuclei and do the RF due to the charges
 
      NU = 0
 
      DO 14 J = 1,NUCIND ! loop only over the QM nuclei (NUCIND). This is the most correct way ... 
        NU = NU + 1      ! The counter could just be over J instead but I keep the refernce to NU
        DO 15 K = 1, ISYTP  ! exclude the QM system here, i.e. begin the loop from 1
          DO 21 LM = NSYSBG(K), NSYSED(K)
            DO 16 L = 1, NCTOT
 
              DIST = D0
              DISM1 = D0
              DISM3 = D0
 
              IF ( (ISUBSY(L) .EQ. LM) .AND. 
     &             (ISUBSI(L) .LE. NSISY(K)) ) THEN

                DIST = (CORD(1,J)-CORD(1,L))**2 +
     &                 (CORD(2,J)-CORD(2,L))**2 +
     &                 (CORD(3,J)-CORD(3,L))**2
                       DIST = SQRT(DIST)
                       DISM1 = 1/DIST
                       DISM3 = DISM1**3
 
                       DO 18 M = 1, 3
                         RFQM(M,NU) = RFQM(M,NU) +
     &                   CHARGE(L)*(CORD(M,J)-CORD(M,L))*DISM3
  18                   CONTINUE
 
              END IF
  16        CONTINUE
  21      CONTINUE
  15    CONTINUE
  14  CONTINUE
 
      IF (NUCIND .NE. NU) THEN
        CALL QUIT('Something is wrong w. the no. of nuclei in RFQMMM')
      ENDIF
 
      WRITE(LUPRI,*) 'The QM coordinates:'
      DO 22 JK=1,NU
        WRITE(LUPRI,*) CORD(1,JK),CORD(2,JK),CORD(3,JK)
  22  CONTINUE
      WRITE(LUPRI,'(/A)') ' The reaction field due to the charges:'
      DO 23 JK=1,NU
        WRITE(LUPRI,*) RFQM(1,JK),RFQM(2,JK),RFQM(3,JK)
  23  CONTINUE

C     Finally, we set the charges to zero if OLDTG
 
      IF (OLDTG) THEN
       DO 876 I= 1, ISYTP
          DO 877 J = NSYSBG(I), NSYSED(I)
            DO 878 K = 1, NCTOT
               IF (ISUBSY(K) .EQ. J) CHARGE(K) = D0
  878       CONTINUE
  877    CONTINUE
  876   CONTINUE
      END IF

C     Now continue with the contribution due to the induced dipole moments

C     1. get induced dipole moments from files:
 
      LUMOM = -1
      CALL QM3_GET3(LUMOM,'QMIND   ',MYIND,NCOMS)
C
      NUM = 0
      NATOM = 0
      DO 30 I = 1, NCTOT
        IF ( (ISUBSY(I) .EQ. 0) .AND.
     &       (ISUBSI(I) .LE. NSISY(0)) ) THEN  ! we have a QM nuclei
           NATOM = NATOM + 1
           NUM = 0
           DO 40 I2 =1, ISYTP ! loop only over MM systens
           IF (MDLWRD(I2)(1:5) .NE. 'SPC_E') GOTO 40
             DO 50 J = NSYSBG(I2), NSYSED(I2)
               DO 60 L=1,NCTOT
                 IF ( (ISUBSY(L) .EQ. J) .AND.
     &                 (ISUBSI(L) .GT. NSISY(I2)) ) THEN
                   TEXX = 0.0D0
                   TEXY = 0.0D0
                   TEXZ = 0.0D0
                   TEYX = 0.0D0
                   TEYY = 0.0D0
                   TEYZ = 0.0D0
                   TEZX = 0.0D0
                   TEZY = 0.0D0
                   TEZZ = 0.0D0
                   NUM = NUM +1
                   XCOOR = CORD(1,I) - CORD(1,L)
                   YCOOR = CORD(2,I) - CORD(2,L)
                   ZCOOR = CORD(3,I) - CORD(3,L)
                   R2 = XCOOR*XCOOR + YCOOR*YCOOR + ZCOOR*ZCOOR
                   R5 = R2*R2*SQRT(R2)
                   TEXX = (D3*XCOOR*XCOOR-R2)/R5
                   TEXY = (D3*XCOOR*YCOOR)/R5
                   TEXZ = (D3*XCOOR*ZCOOR)/R5
                   TEYX = (D3*YCOOR*XCOOR)/R5 
                   TEYY = (D3*YCOOR*YCOOR-R2)/R5
                   TEYZ = (D3*YCOOR*ZCOOR)/R5
                   TEZX = (D3*ZCOOR*XCOOR)/R5
                   TEZY = (D3*ZCOOR*YCOOR)/R5
                   TEZZ = (D3*ZCOOR*ZCOOR-R2)/R5

                   RFQM(1,NATOM) = RFQM(1,NATOM) + MYIND(1,NUM)*TEXX 
     &                                           + MYIND(2,NUM)*TEXY
     &                                           + MYIND(3,NUM)*TEXZ

                   RFQM(2,NATOM) = RFQM(2,NATOM) + MYIND(1,NUM)*TEYX
     &                                           + MYIND(2,NUM)*TEYY
     &                                           + MYIND(3,NUM)*TEYZ

                   RFQM(3,NATOM) = RFQM(3,NATOM) + MYIND(1,NUM)*TEZX
     &                                           + MYIND(2,NUM)*TEZY
     &                                           + MYIND(3,NUM)*TEZZ
                 END IF
 60            CONTINUE
 50          CONTINUE
 40        CONTINUE
        END IF
 30   CONTINUE

      IF (NUCIND .NE. NATOM) THEN
        CALL QUIT('Something is wrong w. the no. of nuclei in RFQMMM')
      ENDIF

      WRITE(LUPRI,'(/A)') ' The QM coordinates:'
      DO 24 JK=1,NATOM
        WRITE(LUPRI,*) CORD(1,JK),CORD(2,JK),CORD(3,JK)
  24  CONTINUE
      WRITE(LUPRI,'(/A)') ' The total reaction field:'
      DO 25 JK=1,NATOM
        WRITE(LUPRI,*) RFQM(1,JK),RFQM(2,JK),RFQM(3,JK)
  25  CONTINUE

      DO 26 JK=1,NATOM
        DO IPRPC = 1,3
           IF (IPRPC .EQ. 1) LAPRPC = 'XREFIELD'
           IF (IPRPC .EQ. 2) LAPRPC = 'YREFIELD'
           IF (IPRPC .EQ. 3) LAPRPC = 'ZREFIELD'
           CALL WRIPRO(RFQM(IPRPC,JK),"SCF/DFT   ",1,
     *                 LAPRPC,LAPRPC,LAPRPC,LAPRPC,
     *                 0.0D0,0.0D0,0.0D0,1,0,0,0)
        END DO
  26  CONTINUE
      END
C***************************************************************************************************************************
C  /* Deck pcmfield */
      SUBROUTINE PCMFIELD(EFPCM)
C
C     Calculates the electric field from the PCM charges at the
C     polarizable points
C
#include "implicit.h"
#include "priunit.h"
#include "qm3.h"
#include "mxcent.h"
#include "inftap.h"
#include "pcmdef.h"
#include "orgcom.h"
#include "pcm.h"
#include "pcmlog.h"
#include "nuclei.h"
 
      PARAMETER  (D0 = 0.0D0)
      REAL*8     EFPCM(3,MXQM3)

      DO 19 I = 1, MXQM3
         DO 20 J = 1, 3
            EFPCM(J,I) = D0
   20    CONTINUE
   19 CONTINUE



      LM = 0
      L = NUSITE + NUALIS(0)
      ITEST = 0

      DO 520 I = 1, ISYTP
        IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
          DO 521 J = NSYSBG(I), NSYSED(I)
            DO 522 K = 1, NUALIS(I)
              LM = LM + 1

              IATNOW = NUCIND + L + LM

C             This is a MM polarazable site
              XCORD = CORD(1,IATNOW)
              YCORD = CORD(2,IATNOW)
              ZCORD = CORD(3,IATNOW)

C             Loop over the contribution from the PCM tesseras
              DO 523 M = 1, NTS
                PCMCHARGE = QSN(M)+QSE(M)+QSMM(M)
                DIST2 = (XCORD-XTSCOR(M))**2
     *               + (YCORD-YTSCOR(M))**2
     *               + (ZCORD-ZTSCOR(M))**2
                DIST = SQRT(DIST2)
                DISTM1 = 1.0/DIST
                DISTM3 = DISTM1**3
                EFPCM(1,LM) = EFPCM(1,LM) 
     *                      + PCMCHARGE*(XCORD-XTSCOR(M))*DISTM3
                EFPCM(2,LM) = EFPCM(2,LM) 
     *                      + PCMCHARGE*(YCORD-YTSCOR(M))*DISTM3
                EFPCM(3,LM) = EFPCM(3,LM) 
     *                      + PCMCHARGE*(ZCORD-ZTSCOR(M))*DISTM3
  523         CONTINUE
  522       CONTINUE
  521     CONTINUE
        END IF
  520 CONTINUE

      IF (LM.NE.NCOMS) THEN
        CALL QUIT('Wrong number of MM pol. sites in PCMFIELD')
      END IF

      END
C*********************************************************************************
C  /* Deck GIVEMY */
      SUBROUTINE GIVEMY(MYIND)
C
C     Transfer the induced dipole moments to pcm
C
#include "implicit.h"
#include "priunit.h"
#include "qm3.h"
#include "mxcent.h"
#include "inftap.h"
#include "pcmdef.h"
#include "orgcom.h"
#include "pcm.h"
#include "pcmlog.h"
#include "nuclei.h"
 
      PARAMETER  (D0 = 0.0D0)
      REAL*8 MYIND(3,MXQM3)

      LM = 0
      L = NUSITE + NUALIS(0)

      DO 520 I = 1, ISYTP
        IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
          DO 521 J = NSYSBG(I), NSYSED(I)
            DO 522 K = 1, NUALIS(I)
              LM = LM + 1

              IATNOW = NUCIND + L + LM

C             This is a MM polarazable site
              XMMMY(LM) = CORD(1,IATNOW)
              YMMMY(LM) = CORD(2,IATNOW)
              ZMMMY(LM) = CORD(3,IATNOW)
              MMMYX(LM) = MYIND(1,LM)
              MMMYY(LM) = MYIND(2,LM)
              MMMYZ(LM) = MYIND(3,LM)
             IF (L.GT.MXQ) THEN
               CALL QUIT('Too many pols.in MM. Increase MXQ')
             ENDIF
  522       CONTINUE
  521     CONTINUE
        END IF
  520 CONTINUE

      IF (LM.NE.NCOMS) THEN
        CALL QUIT('Wrong number of MM pol. sites in PCMFIELD')
      END IF

      END

C**************************************************************************************
C  /* Deck MMPCMINIT */
      SUBROUTINE MMPCMINIT()
C
C     Initializes ICQM3 array keeping information on the convergence of
C     the response vectors
C
C     1: initialized 
C     0: not converged
C    -1: converged 

#include "implicit.h"
#include "qm3.h"

      DO 100 I=1,NSTATES
         ICQM3(I) = 1
  100 CONTINUE 

      END
C***************************************************************************************

C*************************************************************
C  /* Deck MMPCM_FIXMOM */
      SUBROUTINE MMPCM_FIXMOM()
C*************************************************************
C
#include "implicit.h"
#include "priunit.h"
#include "qm3.h"
#include "mxcent.h"
#include "nuclei.h"
#include "inftap.h"
C
C
      REAL*8 MYIND(3,MXQM3)
C
      LUMOM = -1
      CALL QM3_GET3(LUMOM,'QMIND   ',MYIND,NCOMS)
      CALL GIVEMY(MYIND)
C
      END
C-------------------------------------------------------------------------
C  /* Deck qmmmin */
      SUBROUTINE QMMMIN(WORD)
C-------------------------------------------------------------------------
C
C INPUT FROM ABACUS TO QMMM CALCULATIONS 
C
#include "implicit.h"
#include "priunit.h"
#include "gnrinf.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "mxsymm.h"
#include "qm3.h"
#include "qmmm.h"
#include "codata.h"
C
      PARAMETER ( NTABLE = 26, D1 = 1.0D0)
      LOGICAL NEWDEF,RCUT
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER TSTCHA*1,WORD2*6,WORD3*6
      CHARACTER*80 INLINE
      CHARACTER*2 UNITS
#include "cbirea.h"
#include "ccom.h"
C
!                  1          2          3          4
      DATA TABLE /'.QMMM  ', '.SPLDIP', '.FIXDIP', '.MMDAMP',
!                  5          6          7          8
     *            '.MMPROP', '.SPLNMR', '.SELPOL', '.QMDAMP',
!                  9          10         11         12
     *            '.LGSPOL', '.MMMAT ', '.MMITER', '.ITERIN',
!                  13         14         15         16
     *            '.MMDIIS', '.LCLOSE', '.RCUTMM', '.FDSTEP',
!                  17         18         19         20
     *            '.ISOPOL', '.PRINT ', '.MMREST', '.RELMAT',
!                  21         22         23         24
     *            '.ZERODI', '.ZEROQU', '.NOPOL ', '.NOMB  ', 
!                  25         26
     *            '.MMTIME', '.NEWEXC'/
 
      CALL QENTER('QMMMIN')
 
C     Initializing QMMM keywords
 
      QMMM     = .FALSE.   ! Do a NewType QMMM calculation
      SPLDIP   = .FALSE.   ! Split ind. dips. into components from electrons, permanent MM mult. and ind. multipoles.
      FIXDIP   = .FALSE.   ! Do not calc. ind. dipoles read from file. Only for GS.
      MMDAMP   = .FALSE.   ! Damping in relay tensor
      CONMAT   = .TRUE.    ! Indicates if [ALPHA^(-1) - T] matrix is to be constructed.
C                            This is but true here and then after the
C                            first construction this is set to false
      MMPROP   = .FALSE.   ! Calculates properties of the MM environment. 
      SPLNMR   = .FALSE.   ! Calculate London contribution to gradient
C                            seperately from permanent and induced dipoles. For testing.
      QMDAMP   = .FALSE.   ! Include damping of electric field between QM and MM
      LGSPOL   = .FALSE.   ! Skip the response part related to Efield
C      SELPOL   = .FALSE.   ! Skip selected systems response related to the electric field. Response
      MMMAT    = .TRUE.    ! Default is matrix inversion (keyword MATINV)
      MMITER   = .FALSE.   ! Logical for iterative determination of ind. dipoles (keyword ITER)
      MXMMIT   = 100       ! Max. iterations in Jacobi MM induced dipoles solver. 
      THMMIT   = 1.0D-10   ! Threshold in Jacobi/DIIS MM induced dipoles solver (for norm error vector in the induced dipoles))
      MMDIIS   = .FALSE.   ! Use DIIS solver for induced dipoles
      MXMMDI   = 5         ! Maximum vectors to be used in MMDIIS
      LCLOSE   = .FALSE.   ! Logical for investigating close QM-MM distances
      IQMMCC   = -1        ! Option for close QM-MM distance. 
      RQMMMC   = 2.20      ! Close distance QM-MM. Do something if actual distance is less than RQMMMC. In au. 
      RCUTMM   = 1000.0D0  ! Cut-off for QM-MM and MM-MM interactions
      RCUT     = .FALSE.   ! No cut-off
      DELFLD   = 0.01      ! Keywd FDSTEP: FD step in calculation of electric field Hessian integrals
      ZERODI   = .FALSE.   ! Zero the dipole moments
      ZEROQU   = .FALSE.   ! Zero the quadrupole moments
      ISOPOL   = .FALSE.   ! Make anisotropic polarizabilities isotropic
      NOPOL    = .FALSE.   ! Zero the static polarizabilities
      IPQMMM   = 1         ! QMMM print level
      MMREST   = .FALSE.   ! Restart QM/MM calculation (only useful for reusing relay matrix atm.)
      RELMAT   = .FALSE.   ! Construct relay matrix, save to QMMMIM and quit
      NOMB     = .FALSE.   ! Remove many-body effects, i.e. no induced dipole - induced dipole interactions
      MMTIME   = .FALSE.   ! Print wall times used in different parts of the QMMM module
      NEWEXC   = .FALSE.   ! Use new type of exclusion list

      NEWDEF = (WORD .EQ. '*QMMM  ')
      ICHANG = 0
C
      IF (NEWDEF) THEN
        WORD1 = WORD
C
  100   CONTINUE
C
        READ (LUCMD, '(A7)') WORD
        CALL UPCASE(WORD)
        PROMPT = WORD(1:1)
C
        IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
          GO TO 100
        ELSE IF (PROMPT .EQ. '.') THEN
          ICHANG = ICHANG + 1
          DO 200 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
              GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
     &               19,20,21,22,23,24,25,26), I
            END IF
  200     CONTINUE
          IF (WORD .EQ. '.OPTION') THEN
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            GO TO 100
          END IF
C
          WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *                             '" not recognized for *QMMM  .'
          CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
          CALL QUIT('Illegal keyword in REAINP.')
C
    1     CONTINUE
            QMMM = .TRUE.
            NYQMMM = .TRUE.
            LGFILE = .FALSE. !sne assignment. 
            LHFFIL = .FALSE.
          GO TO 100
C
    2     CONTINUE
            SPLDIP = .TRUE.
          GO TO 100
C
    3     CONTINUE
            FIXDIP = .TRUE.
          GO TO 100
C
    4     CONTINUE
            MMDAMP = .TRUE.
          GO TO 100
C
    5     CONTINUE
            MMPROP = .TRUE.
          GO TO 100
C
    6     CONTINUE
            SPLNMR = .TRUE.
          GO TO 100
C
    7     CONTINUE
c            SELPOL = .TRUE.
c            READ(LUCMD,*) 
          GO TO 100
C
    8     CONTINUE
            QMDAMP = .TRUE.
            READ(LUCMD,*) IDAMP
            IF (IDAMP .EQ. 3) THEN
              READ(LUCMD,*) NQMNUC
              DO 999 J = 1,NQMNUC
                READ(LUCMD,*) QMPOL(J)
  999         CONTINUE
            ELSE 
              READ(LUCMD,*) ADAMP
            ENDIF
          GO TO 100

    9     CONTINUE
            LGSPOL = .TRUE.
          GO TO 100

   10     CONTINUE
            MMMAT = .TRUE.
            MMITER = .FALSE.
          GO TO 100

   11     CONTINUE
            MMMAT = .FALSE.
            MMITER = .TRUE.
          GO TO 100

   12     CONTINUE
            READ(LUCMD,*) MXMMIT, THMMIT
          GO TO 100

   13     CONTINUE
            MMDIIS = .TRUE.
            READ(LUCMD,*) MXMMDI
          GO TO 100

   14     CONTINUE
            LCLOSE = .TRUE.
C           See if further input is provided for this option        
            READ (LUCMD,'(A7)') WORD
            CALL UPCASE(WORD)
            BACKSPACE(LUCMD)
            IF (WORD(1:1).NE.'.' .AND. WORD(1:1).NE.'*'
     &                           .AND. WORD(1:1).NE.'!' ) THEN
              READ (LUCMD,*) RQMMMC, UNITS, IQMMMC
              CALL UPCASE(UNITS)
              IF (UNITS .EQ. 'AA') THEN
                RQMMMC = RQMMMC/XTANG
              ELSE IF (UNITS .EQ. 'AU') THEN
                RQMMMC = RQMMMC
              ELSE
                CALL QUIT('Unknown units for LCLOSE specification')
              ENDIF
            END IF
          GO TO 100

   15     CONTINUE
            READ(LUCMD,*) RCUTMM, UNITS 
            CALL UPCASE(UNITS)
            IF (UNITS .EQ. 'AA') THEN
              RCUTMM = RCUTMM/XTANG
            ELSE IF (UNITS .EQ. 'AU') THEN
              RCUTMM = RCUTMM
            ELSE
              CALL QUIT('Unknown units for RCUTMM specification')
            ENDIF
            RCUT = .TRUE.
          GO TO 100

   16     CONTINUE
            READ(LUCMD,*) DELFLD
          GO TO 100

   17     CONTINUE
            ISOPOL = .TRUE.
          GO TO 100

   18     CONTINUE
            READ(LUCMD,*) IPQMMM
          GO TO 100

   19     CONTINUE
            MMREST = .TRUE.
            IF (MMMAT) THEN
              INQUIRE(FILE='QMMMIM',EXIST=CONMAT)
              IF (CONMAT) THEN
                CONMAT = .FALSE.
              ELSE
                WRITE (LUPRI,'(A)')
                WRITE (LUPRI,'(A)') 'File containing the relay matrix'//
     &                              ' not found!'
                WRITE (LUPRI,'(A)')
                WRITE (LUPRI,'(A)') 'Remember to copy it to the'//
     &                              ' scratch directory either'//
     &                              ' manually or with the'
                WRITE (LUPRI,'(A)') 'dalton script using the "-qmmm"'//
     &                              ' option.'
                WRITE (LUPRI,'(A)')
                CALL QUIT('ERROR: QMMMIM file not found')
              ENDIF
            ENDIF
          GO TO 100

   20     CONTINUE
            RELMAT = .TRUE.
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Relay matrix will be constructed and saved'
            WRITE(LUPRI,*) 'in QMMMIM.'
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Skipping everything else!'
          GO TO 100

   21     CONTINUE
            ZERODI = .TRUE.
            ZEROQU = .TRUE.
          GO TO 100

   22     CONTINUE
            ZEROQU = .TRUE.
          GO TO 100

   23     CONTINUE
            NOPOL = .TRUE.
          GO TO 100

   24     CONTINUE
            NOMB = .TRUE.
          GO TO 100

   25     CONTINUE
            MMTIME = .TRUE.
          GO TO 100

   26     CONTINUE
            NEWEXC = .TRUE.
          GO TO 100

        ELSE IF (PROMPT .EQ. '*') THEN
          GO TO 300
        ELSE
          WRITE (LUPRI,'(/3A/)') ' Prompt "',WORD,
     *                             '" not recognized for *QMMM  .'
          CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
          CALL QUIT('Illegal prompt for *QMMM  .')
        END IF
      END IF
C
  300 CONTINUE
C
      IF (ICHANG .GT. 0) THEN
        CALL HEADER('Changes of defaults for *QMMM  :',0)
C
        WRITE(LUPRI,'(1X,A1,30A,A1)') '+',('-',J=1,18),'+'
        WRITE (LUPRI,'(1X,A2,A6,A3,A7,A2)')
     *         '| ',' WORD:',' | ','CHANGE:',' |'
        WRITE(LUPRI,'(1X,A1,30A,A1)') '+',('-',J=1,18),'+'
C
        IF (RELMAT) THEN
                  WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ','RELMAT',' | ',RELMAT,' |'
          GOTO 400
        ENDIF
        IF (QMMM) WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ','  QMMM',' | ',QMMM,' |'
        IF (SPLDIP) WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ','SPLDIP',' | ',SPLDIP,' |'
        IF (FIXDIP) WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ','FIXDIP',' | ',FIXDIP,' |'
        IF (MMDAMP) WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ','MMDAMP',' | ',MMDAMP,' |'
        IF (MMPROP) WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ','MMPROP',' | ',MMPROP,' |'
        IF (SPLNMR) WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ','SPLNMR',' | ',SPLNMR,' |'
        IF (QMDAMP) WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ','QMDAMP',' | ',QMDAMP,' |'
        IF (LGSPOL) WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ','LGSPOL',' | ',LGSPOL,' |'
        IF (MMITER) WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ','MMITER',' | ',MMITER,' |'
        IF (MMDIIS) WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ','MMDIIS',' | ',MMDIIS,' |'
        IF (LCLOSE) WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ','LCLOSE',' | ',LCLOSE,' |'
        IF (RCUT) WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ','RCUT',' | ',RCUT,' |'
        IF (MMREST) WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ','MMREST',' | ',MMREST,' |'
        IF (NOMB) WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ','NOMB',' | ',NOMB,' |'
        IF (MMTIME) WRITE (LUPRI,'(1X,A2,A6,A3,L7,A2)')
     *         '| ','MMTIME',' | ',MMTIME,' |'
  400   WRITE(LUPRI,'(A2,30A1)') ' +',('-',J=1,18),'+'
        WRITE(LUPRI,*)
      END IF

      IF (MMITER) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) 'Induced MM dipoles are solved iteratively'
        WRITE(LUPRI,*) 'Max. number of iterations: ',MXMMIT
        WRITE(LUPRI,*) 'Threshold: ',THMMIT
        IF (MMDIIS) THEN
          WRITE(LUPRI,*) 'DIIS is invoked for induced dipoles'
          WRITE(LUPRI,*) 'Maximum number of DIIS vectors: ',MXMMDI
        ENDIF
        WRITE(LUPRI,*)
      ENDIF

      IF (LCLOSE) THEN
        WRITE(LUPRI,*) 'Close QM-MM distance threshold (au) ', RQMMMC
        IF (IQMMMC .GT. 0) THEN
          WRITE(LUPRI,*) 'Close QM-MM distance criteria: ', IQMMMC
        ENDIF
      ENDIF

      IF ( (SPLDIP) .AND. (.NOT. MMMAT) ) THEN
        WRITE(LUPRI,*) 'Split only allowed if matrix inversion'
        CALL QUIT('Error in QMMM input')
      ENDIF

      IF (DELFLD .NE. 0.01) THEN
        WRITE(LUPRI,*) 'Step-length in calc. of Field-Hessian integrals'
        WRITE(LUPRI,*) DELFLD
      ENDIF

      IF (ZERODI) THEN
        WRITE(LUPRI,*) 'Dipoles are reset to zero'
      ENDIF 

      IF (ZEROQU) THEN
        WRITE(LUPRI,*) 'Quadrupoles are reset to zero'
      ENDIF 

      IF (ISOPOL) THEN
        WRITE(LUPRI,*) 'Anisotropic polarizabilites are made isotropic'
      ENDIF 

      IF (NOPOL) THEN
        WRITE(LUPRI,*) 'Static polarizabilites are reset to zero'
      ENDIF

      IF (NOMB) THEN
        WRITE(LUPRI,*) 'No induced dipole - induced dipole interactions'
      ENDIF

      CALL QEXIT('QMMMIN')
      RETURN
      END
C******************************************************************************
C  /* Deck readmm */
      SUBROUTINE READMM(WORK,LWORK)
C------------------------------------------------------------------------
C
C --- General Input for QMMM electrostatic potentials---
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "mxcent.h"
#include "abainf.h"
#include "inftap.h"
#include "nuclei.h"
#include "qmmm.h"
#include "codata.h"
C
      CHARACTER*180 KWord
      CHARACTER*2 UNITS
      CHARACTER*25 POLWRD
      DIMENSION WORK(LWORK)
      LOGICAL LFIND,LHIT,LREAD

      CALL QENTER('READMM')

      LREAD = .FALSE.

C     First we zero everything
      DO 099 J = 1,MMCENT
        MMCORD(1,J) = 0.0D0
        MMCORD(2,J) = 0.0D0
        MMCORD(3,J) = 0.0D0
        MUL0MM(J)   = 0.0D0
        MUL1MM(1,J) = 0.0D0
        MUL1MM(2,J) = 0.0D0
        MUL1MM(3,J) = 0.0D0
        MUL2MM(1,J) = 0.0D0
        MUL2MM(2,J) = 0.0D0
        MUL2MM(3,J) = 0.0D0
        MUL2MM(4,J) = 0.0D0
        MUL2MM(5,J) = 0.0D0
        MUL2MM(6,J) = 0.0D0
        POLMM(1,J)  = 0.0D0
        POLMM(2,J)  = 0.0D0
        POLMM(3,J)  = 0.0D0
        POLMM(4,J)  = 0.0D0
        POLMM(5,J)  = 0.0D0
        POLMM(6,J)  = 0.0D0
        POLIMM(J)   = 0.0D0
        LMMNUL(J)   = .FALSE.
  099 CONTINUE

      LUQMMM=-1
      CALL GPOPEN(LUQMMM,'POTENTIAL.INP','OLD',' ',
     &            'FORMATTED',IDUMMY,.FALSE.)
      REWIND(LUQMMM)
      READ(LUQMMM,*) UNITS

      CALL UPCASE(UNITS)

      IF (UNITS .EQ. 'AA') THEN
        SCAL = XTANG
      ELSE IF (UNITS .EQ. 'AU') THEN
        SCAL = 1.0D0
      ELSE
         CALL QUIT('Unknown units in POTENTIAL.INP')
      ENDIF

      Read(LUQMMM,'(A180)') KWord 
      KWord(168:180)='-2 -2 -2 -2'
      Read(KWord,*)MMCENT,NMULT,IPOLTP,NEXLST,NELEME

      IF (MMCENT .GT. MXMMCT) THEN 
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) 'MXMMCT:', MXMMCT
        WRITE(LUPRI,*) 'No. of centers read from qmmm.inp:', MMCENT
        CALL QUIT('Increase MXMMCT')
      ENDIF

      IF (NEXLST .GT. MXEXCL) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) 'MXEXCL:', MXEXCL
        WRITE(LUPRI,*) 'Length of excl. list from qmmm.inp:', NEXLST
        CALL QUIT('Increase NEXLST')
      ENDIF

C     default setting: 
C     - order of multipoles = 0 (charges)    
C     - polarization: No
C     - length of excluded list = 0

      IF (NMULT .EQ. -2) NMULT = 0
      IF (IPOLTP .EQ. -2) IPOLTP = 0
      IF (NEXLST .EQ. -2) NEXLST = 0
      IF (NELEME .EQ. -2) NELEME = 0

      IF (NOPOL) IPOLTP = 0
      IF (IPOLTP .EQ. 0) POLWRD = ' No polarization         '
      IF (IPOLTP .EQ. 1) POLWRD = ' Isotropic polarization  '
      IF (IPOLTP .EQ. 2) POLWRD = ' Anisotropic polarization'

      WRITE(LUPRI,*)
      WRITE(LUPRI,*) ' *********************************** '
      WRITE(LUPRI,*) ' QMMM electrostatic potential: '
      WRITE(LUPRI,*) ' Multipole order               ', NMULT
      WRITE(LUPRI,*) POLWRD
      WRITE(LUPRI,*) ' *********************************** '

      IF (NMULT. GT. 2) THEN
        CALL QUIT('Multipole order only impl. for 0,1,2')
      ENDIF

      IF ((IPOLTP. NE. 0) .AND. (IPOLTP. NE. 1)
     *                    .AND. (IPOLTP. NE. 2)) THEN
        CALL QUIT('Polarizability order only impl. for 0,1,2')
      ENDIF

C     Only charges ...  No polarizabilities 

      IF ( (NMULT .EQ. 0) .AND. (IPOLTP .EQ. 0) ) THEN
        LREAD = .TRUE.
        IF (NELEME .EQ. 0) THEN
          DO 100 J = 1,MMCENT
            READ(LUQMMM,*) (EXLIST(I,J), I=1,NEXLST),
     *                     MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                     MUL0MM(J)
            MMCORD(1,J) = MMCORD(1,J)/SCAL
            MMCORD(2,J) = MMCORD(2,J)/SCAL
            MMCORD(3,J) = MMCORD(3,J)/SCAL
  100     CONTINUE
        ELSE
          DO 300 J = 1,MMCENT
            READ(LUQMMM,*) (EXLIST(I,J), I=1,NEXLST),ELEME(J),
     *                     MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                     MUL0MM(J)
            MMCORD(1,J) = MMCORD(1,J)/SCAL
            MMCORD(2,J) = MMCORD(2,J)/SCAL
            MMCORD(3,J) = MMCORD(3,J)/SCAL
  300     CONTINUE
        ENDIF
      ENDIF

C     Charges + dipoles ...  No polarizabilities   

      IF ( (NMULT .EQ. 1) .AND. (IPOLTP .EQ. 0) ) THEN
        LREAD = .TRUE.
        IF (NELEME .EQ. 0) THEN 
          DO 101 J = 1,MMCENT
            READ(LUQMMM,*) (EXLIST(I,J), I=1,NEXLST),
     *                     MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                     MUL0MM(J),
     *                     MUL1MM(1,J),MUL1MM(2,J),MUL1MM(3,J)
            MMCORD(1,J) = MMCORD(1,J)/SCAL
            MMCORD(2,J) = MMCORD(2,J)/SCAL
            MMCORD(3,J) = MMCORD(3,J)/SCAL
  101     CONTINUE
        ELSE
          DO 301 J = 1,MMCENT
            READ(LUQMMM,*) (EXLIST(I,J), I=1,NEXLST),ELEME(J),
     *                     MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                     MUL0MM(J),
     *                     MUL1MM(1,J),MUL1MM(2,J),MUL1MM(3,J)
            MMCORD(1,J) = MMCORD(1,J)/SCAL
            MMCORD(2,J) = MMCORD(2,J)/SCAL
            MMCORD(3,J) = MMCORD(3,J)/SCAL
  301     CONTINUE
        ENDIF
      ENDIF

C     Charges + dipoles + quadrupoles ...  No polarizabilities   

      IF ( (NMULT .EQ. 2) .AND. (IPOLTP .EQ. 0) ) THEN
        LREAD = .TRUE.
        IF (NELEME .EQ. 0) THEN
          DO 102 J = 1,MMCENT
            READ(LUQMMM,*) (EXLIST(I,J), I=1,NEXLST),
     *                     MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                     MUL0MM(J),
     *                     MUL1MM(1,J),MUL1MM(2,J),MUL1MM(3,J),
     *                     MUL2MM(1,J),MUL2MM(2,J),MUL2MM(3,J),
     *                     MUL2MM(4,J),MUL2MM(5,J),MUL2MM(6,J)
            MMCORD(1,J) = MMCORD(1,J)/SCAL
            MMCORD(2,J) = MMCORD(2,J)/SCAL
            MMCORD(3,J) = MMCORD(3,J)/SCAL
  102     CONTINUE
        ELSE
          DO 302 J = 1,MMCENT
            READ(LUQMMM,*) (EXLIST(I,J), I=1,NEXLST),ELEME(J),
     *                     MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                     MUL0MM(J),
     *                     MUL1MM(1,J),MUL1MM(2,J),MUL1MM(3,J),
     *                     MUL2MM(1,J),MUL2MM(2,J),MUL2MM(3,J),
     *                     MUL2MM(4,J),MUL2MM(5,J),MUL2MM(6,J)
            MMCORD(1,J) = MMCORD(1,J)/SCAL
            MMCORD(2,J) = MMCORD(2,J)/SCAL
            MMCORD(3,J) = MMCORD(3,J)/SCAL
  302     CONTINUE
        ENDIF
      ENDIF

C     Only charges + anisotropic polarizabiblities 

      IF ( (NMULT .EQ. 0) .AND. (IPOLTP .EQ. 2) ) THEN
        LREAD = .TRUE.
        IF (NELEME .EQ. 0) THEN
          DO 103 J = 1,MMCENT
            READ(LUQMMM,*) (EXLIST(I,J), I=1,NEXLST),
     *                     MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                     MUL0MM(J),
     *                     POLMM(1,J),POLMM(2,J),POLMM(3,J),
     *                     POLMM(4,J),POLMM(5,J),POLMM(6,J)
            MMCORD(1,J) = MMCORD(1,J)/SCAL
            MMCORD(2,J) = MMCORD(2,J)/SCAL
            MMCORD(3,J) = MMCORD(3,J)/SCAL
  103     CONTINUE
       ELSE
          DO 303 J = 1,MMCENT
            READ(LUQMMM,*) (EXLIST(I,J), I=1,NEXLST),ELEME(J),
     *                     MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                     MUL0MM(J),
     *                     POLMM(1,J),POLMM(2,J),POLMM(3,J),
     *                     POLMM(4,J),POLMM(5,J),POLMM(6,J)
            MMCORD(1,J) = MMCORD(1,J)/SCAL
            MMCORD(2,J) = MMCORD(2,J)/SCAL
            MMCORD(3,J) = MMCORD(3,J)/SCAL
  303     CONTINUE
       ENDIF
      ENDIF

C     Charges + dipoles  + anisotropic polarizabiblities   

      IF ( (NMULT .EQ. 1) .AND. (IPOLTP .EQ. 2) ) THEN
        LREAD = .TRUE.
        IF (NELEME .EQ. 0) THEN
          DO 104 J = 1,MMCENT
            READ(LUQMMM,*) (EXLIST(I,J), I=1,NEXLST),
     *                     MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                     MUL0MM(J),
     *                     MUL1MM(1,J),MUL1MM(2,J),MUL1MM(3,J),
     *                     POLMM(1,J),POLMM(2,J),POLMM(3,J),
     *                     POLMM(4,J),POLMM(5,J),POLMM(6,J)
            MMCORD(1,J) = MMCORD(1,J)/SCAL
            MMCORD(2,J) = MMCORD(2,J)/SCAL
            MMCORD(3,J) = MMCORD(3,J)/SCAL
  104     CONTINUE
        ELSE
          DO 304 J = 1,MMCENT
            READ(LUQMMM,*) (EXLIST(I,J), I=1,NEXLST),ELEME(J),
     *                     MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                     MUL0MM(J),
     *                     MUL1MM(1,J),MUL1MM(2,J),MUL1MM(3,J),
     *                     POLMM(1,J),POLMM(2,J),POLMM(3,J),
     *                     POLMM(4,J),POLMM(5,J),POLMM(6,J)
            MMCORD(1,J) = MMCORD(1,J)/SCAL
            MMCORD(2,J) = MMCORD(2,J)/SCAL
            MMCORD(3,J) = MMCORD(3,J)/SCAL
  304     CONTINUE
        ENDIF
      ENDIF

C     Charges + dipoles + quadrupoles + anisotropic polarizabiblities  

      IF ( (NMULT .EQ. 2) .AND. (IPOLTP .EQ. 2) ) THEN
        LREAD = .TRUE.
        IF (NELEME .EQ. 0) THEN
          DO 105 J = 1,MMCENT
            READ(LUQMMM,*) (EXLIST(I,J), I=1,NEXLST),
     *                     MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                     MUL0MM(J),
     *                     MUL1MM(1,J),MUL1MM(2,J),MUL1MM(3,J),
     *                     MUL2MM(1,J),MUL2MM(2,J),MUL2MM(3,J),
     *                     MUL2MM(4,J),MUL2MM(5,J),MUL2MM(6,J),
     *                     POLMM(1,J),POLMM(2,J),POLMM(3,J),
     *                     POLMM(4,J),POLMM(5,J),POLMM(6,J)
            MMCORD(1,J) = MMCORD(1,J)/SCAL
            MMCORD(2,J) = MMCORD(2,J)/SCAL
            MMCORD(3,J) = MMCORD(3,J)/SCAL
  105     CONTINUE
        ELSE
          DO 305 J = 1,MMCENT
            READ(LUQMMM,*) (EXLIST(I,J), I=1,NEXLST),ELEME(J),
     *                     MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                     MUL0MM(J),
     *                     MUL1MM(1,J),MUL1MM(2,J),MUL1MM(3,J),
     *                     MUL2MM(1,J),MUL2MM(2,J),MUL2MM(3,J),
     *                     MUL2MM(4,J),MUL2MM(5,J),MUL2MM(6,J),
     *                     POLMM(1,J),POLMM(2,J),POLMM(3,J),
     *                     POLMM(4,J),POLMM(5,J),POLMM(6,J)
            MMCORD(1,J) = MMCORD(1,J)/SCAL
            MMCORD(2,J) = MMCORD(2,J)/SCAL
            MMCORD(3,J) = MMCORD(3,J)/SCAL
  305     CONTINUE
        ENDIF
      ENDIF

C     Only charges + isotropic polarizabiblities 

      IF ( (NMULT .EQ. 0) .AND. (IPOLTP .EQ. 1) ) THEN
        LREAD = .TRUE.
        IF (NELEME .EQ. 0) THEN
          DO 106 J = 1,MMCENT
            READ(LUQMMM,*) (EXLIST(I,J), I=1,NEXLST),
     *                     MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                     MUL0MM(J),
     *                     POLIMM(J)
            MMCORD(1,J) = MMCORD(1,J)/SCAL
            MMCORD(2,J) = MMCORD(2,J)/SCAL
            MMCORD(3,J) = MMCORD(3,J)/SCAL
  106     CONTINUE
        ELSE
          DO 306 J = 1,MMCENT
            READ(LUQMMM,*) (EXLIST(I,J), I=1,NEXLST),ELEME(J),
     *                     MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                     MUL0MM(J),
     *                     POLIMM(J)
            MMCORD(1,J) = MMCORD(1,J)/SCAL
            MMCORD(2,J) = MMCORD(2,J)/SCAL
            MMCORD(3,J) = MMCORD(3,J)/SCAL
  306     CONTINUE
        ENDIF
      ENDIF

C     Charges + dipoles  + isotropic polarizabiblities   

      IF ( (NMULT .EQ. 1) .AND. (IPOLTP .EQ. 1) ) THEN
        LREAD = .TRUE.
        IF (NELEME .EQ. 0) THEN
          DO 107 J = 1,MMCENT
            READ(LUQMMM,*) (EXLIST(I,J), I=1,NEXLST),
     *                     MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                     MUL0MM(J),
     *                     MUL1MM(1,J),MUL1MM(2,J),MUL1MM(3,J),
     *                     POLIMM(J)
            MMCORD(1,J) = MMCORD(1,J)/SCAL
            MMCORD(2,J) = MMCORD(2,J)/SCAL
            MMCORD(3,J) = MMCORD(3,J)/SCAL
  107     CONTINUE
        ELSE
          DO 307 J = 1,MMCENT
            READ(LUQMMM,*) (EXLIST(I,J), I=1,NEXLST),ELEME(J),
     *                     MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                     MUL0MM(J),
     *                     MUL1MM(1,J),MUL1MM(2,J),MUL1MM(3,J),
     *                     POLIMM(J)
            MMCORD(1,J) = MMCORD(1,J)/SCAL
            MMCORD(2,J) = MMCORD(2,J)/SCAL
            MMCORD(3,J) = MMCORD(3,J)/SCAL
  307     CONTINUE
        ENDIF
      ENDIF

C     Charges + dipoles + quadrupoles + isotropic polarizabiblities  

      IF ( (NMULT .EQ. 2) .AND. (IPOLTP .EQ. 1) ) THEN
        LREAD = .TRUE.
        IF (NELEME .EQ. 0) THEN
          DO 108 J = 1,MMCENT
            READ(LUQMMM,*) (EXLIST(I,J), I=1,NEXLST), 
     *                     MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                     MUL0MM(J),
     *                     MUL1MM(1,J),MUL1MM(2,J),MUL1MM(3,J),
     *                     MUL2MM(1,J),MUL2MM(2,J),MUL2MM(3,J),
     *                     MUL2MM(4,J),MUL2MM(5,J),MUL2MM(6,J),
     *                     POLIMM(J)
            MMCORD(1,J) = MMCORD(1,J)/SCAL
            MMCORD(2,J) = MMCORD(2,J)/SCAL
            MMCORD(3,J) = MMCORD(3,J)/SCAL
  108     CONTINUE
        ELSE
          DO 308 J = 1,MMCENT
            READ(LUQMMM,*) (EXLIST(I,J), I=1,NEXLST),ELEME(J),
     *                     MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                     MUL0MM(J),
     *                     MUL1MM(1,J),MUL1MM(2,J),MUL1MM(3,J),
     *                     MUL2MM(1,J),MUL2MM(2,J),MUL2MM(3,J),
     *                     MUL2MM(4,J),MUL2MM(5,J),MUL2MM(6,J),
     *                     POLIMM(J)
            MMCORD(1,J) = MMCORD(1,J)/SCAL
            MMCORD(2,J) = MMCORD(2,J)/SCAL
            MMCORD(3,J) = MMCORD(3,J)/SCAL
  308     CONTINUE
        ENDIF
      ENDIF

      CALL GPCLOSE(LUQMMM,'KEEP')

      IF (.NOT. LREAD) THEN
        CALL QUIT('No read option for this potential')
      ENDIF

C     Remove quadrupoles
      IF ( (ZEROQU) .AND. (NMULT .GE. 2) ) THEN
         DO 322 J = 1,MMCENT
            DO 323 I=1,6
               MUL2MM(I,J) = 0.0D0
 323        CONTINUE
 322     CONTINUE 
         NMULT = 1
      ENDIF
C     Remove dipoles
      IF ( (ZERODI) .AND. (NMULT .GE. 1) ) THEN
         DO 320 J = 1,MMCENT
            DO 321 I=1,3
               MUL1MM(I,J) = 0.0D0
 321        CONTINUE
 320     CONTINUE 
         NMULT = 0
      ENDIF
C     Make anisotropic polarizabilities isotropic
      IF (ISOPOL .AND. (IPOLTP .EQ. 2)) THEN
         DO 326 J = 1,MMCENT
            POLIMM(J) = POLMM(1,J)+POLMM(4,J)+POLMM(6,J)
            POLIMM(J) = POLIMM(J)/3.0D0
 326     CONTINUE 
         IPOLTP = 1
      ENDIF

      IOPTMM = 1 
      CALL OUTMM(IOPTMM)

C     Define MM bonds
      NMAX = 0
      IF (NELEME .NE. 0) THEN

        DO 10 J= 1,MMCENT
          IF (ELEME(J) .EQ. 0) GOTO 10
          NBONDS = 0
          RADJ = RADIUS(ELEME(J))
          IF (RADJ .LT. 0.0D0 .AND. ELEME(J) .GT. 0) THEN
              WRITE(LUPRI,*)
     &        'RADIUS FOR MM ATOM WITH ATOMIC NUMBER ',
     &        ELEME(J),' IS UNAVAILABLE, USING 2.0 AA'
              RADJ = 2.0D0
          END IF
          DO 20 I= 1,MMCENT
            IF (ELEME(I) .EQ. 0) GOTO 20
            RADI = RADIUS(ELEME(I))
            IF (RADI .LT. 0.0D0 .AND. ELEME(I) .GT. 0) THEN
              RADI = 2.0D0
            END IF
            IF (I .NE. J) THEN
              DIST2 = (MMCORD(1,J)-MMCORD(1,I))**2 + 
     &                (MMCORD(2,J)-MMCORD(2,I))**2 +
     &                (MMCORD(3,J)-MMCORD(3,I))**2            
              DIST = SQRT(DIST2)*XTANG
              IF (DIST .LT. (1.2D0 * (RADI + RADJ))) THEN
                NBONDS = NBONDS + 1
              END IF
            ENDIF
 20       CONTINUE
          IF (NBONDS .GT. NMAX) NMAX = NBONDS
 10     CONTINUE

        KSTART   = 1
        KEND     = KSTART + (NMAX+1)*MMCENT
        LWORK0   = LWORK  - KEND
        IF (LWORK0 .LT. 0) CALL ERRWRK('READMM',-KEND,LWORK)
        CALL DZERO(WORK(KSTART),NMAX*MMCENT)

        DO 30 J= 1,MMCENT
          IF (ELEME(J) .EQ. 0) GOTO 30
          RADJ = RADIUS(ELEME(J))
          IF (RADJ .LT. 0.0D0 .AND. ELEME(J) .GT. 0) THEN
              WRITE(LUPRI,*)
     &        'RADIUS FOR MM ATOM WITH ATOMIC NUMBER ',
     &        ELEME(J),' IS UNAVAILABLE, USING 2.0 AA'
              RADJ = 2.0D0
          END IF

          IOFF = (J-1)*(NMAX+1)
          WORK(KSTART+IOFF) = J

          LM = 0
          DO 40 I= 1,MMCENT
            IF (ELEME(I) .EQ. 0) GOTO 40
            RADI = RADIUS(ELEME(I))
            IF (RADI .LT. 0.0D0 .AND. ELEME(I) .GT. 0) THEN
              RADI = 2.0D0
            END IF
            IF (I .NE. J) THEN
              DIST2 = (MMCORD(1,J)-MMCORD(1,I))**2 +
     &                (MMCORD(2,J)-MMCORD(2,I))**2 +
     &                (MMCORD(3,J)-MMCORD(3,I))**2
              DIST = SQRT(DIST2)*XTANG
              IF (DIST .LT. (1.2D0 * (RADI + RADJ))) THEN
                LM = LM+1
                WORK(KSTART+IOFF+LM) = I
              END IF
            ENDIF
 40       CONTINUE
 30     CONTINUE

        WRITE(LUPRI,*)
        WRITE(LUPRI,*) '*** MM Topology ***'
        DO 50 I= 1,MMCENT
          WRITE(LUPRI,*) 
     *   (INT(WORK(KSTART+(I-1)*(NMAX+1)+(J-1))),J=1,NMAX+1)
 50     CONTINUE
 
      ENDIF !(NELEME .NE. 0)


C     Do some testing on the QM to MM distances

      IF (LCLOSE .AND. (IQMMMC .EQ. 1)) THEN ! Simple case: put MM elec. parameters to zero

        WRITE(LUPRI,*)
        TEMPQ = 0.0D0
        DO 206 I = 1,NUCIND
          TEMPQ1 = 0.0D0
          DO 207 J = 1,MMCENT
            DIST2 = (CORD(1,I)-MMCORD(1,J))**2 +
     *              (CORD(2,I)-MMCORD(2,J))**2 +
     *              (CORD(3,I)-MMCORD(3,J))**2
            DIST = SQRT(DIST2)
            IF (DIST .LE. RQMMMC) THEN
              WRITE(LUPRI,*) 'QM-MM distance for pair',I,' ',J,' ',DIST,
     *                     ' is small'
              IF (NMULT .GE. 0) THEN
                TEMPQ = TEMPQ + MUL0MM(J)
                MUL0MM(J) = 0.0D0
                LMMNUL(J) = .TRUE.
              ENDIF
              IF (NMULT .GE. 1) THEN
                MUL1MM(1,J) = 0.0D0
                MUL1MM(2,J) = 0.0D0
                MUL1MM(3,J) = 0.0D0
                LMMNUL(J) = .TRUE.
              ENDIF
              IF (NMULT .GE. 2) THEN
                MUL2MM(1,J) = 0.0D0
                MUL2MM(2,J) = 0.0D0
                MUL2MM(3,J) = 0.0D0
                MUL2MM(4,J) = 0.0D0
                MUL2MM(5,J) = 0.0D0
                MUL2MM(6,J) = 0.0D0
                LMMNUL(J) = .TRUE.
              ENDIF
              IF (IPOLTP .EQ. 1) THEN
                POLIMM(J) = 0.0D0
                LMMNUL(J) = .TRUE.
              ELSE IF (IPOLTP .EQ. 2) THEN
                POLMM(1,J) = 0.0D0
                POLMM(2,J) = 0.0D0
                POLMM(3,J) = 0.0D0
                POLMM(4,J) = 0.0D0
                POLMM(5,J) = 0.0D0
                POLMM(6,J) = 0.0D0
                LMMNUL(J) = .TRUE.
              ENDIF
            ENDIF
  207     CONTINUE
  206   CONTINUE

      ENDIF

      IF (LCLOSE .AND. (IQMMMC .EQ. 2)) THEN !Move MM props to nearest other MM site

        TEMPQ = 0.0D0

        KIQM   = 1
        KIMM   = KIQM   + MXMMCT ! we use here just MXMMCT for the length of these temp. arrays
        KEND   = KIMM   + MXMMCT
        LWORK1 = LWORK  - KEND

        IF (LWORK1 .LT. 0) CALL ERRWRK('READMM',-KEND,LWORK)

        CALL DZERO(WORK(KIQM),MXMMCT)
        CALL DZERO(WORK(KIMM),MXMMCT)

C       First find the QM-MM close centers and save those
        LM=0
        DO 208 I = 1,NUCIND
          DO 209 J = 1,MMCENT
            DIST2 = (CORD(1,I)-MMCORD(1,J))**2 +
     *              (CORD(2,I)-MMCORD(2,J))**2 +
     *              (CORD(3,I)-MMCORD(3,J))**2
            DIST = SQRT(DIST2)
            IF (DIST .LE. RQMMMC) THEN
              WRITE(LUPRI,*)
              WRITE(LUPRI,*) 'QM-MM distance for pair',I,' ',J,' ',DIST,
     *                     ' is small'
              WORK(KIQM+LM) = I
              WORK(KIMM+LM) = J
              LM=LM+1  ! actual length of KIMM/KIQM are now LM
              IF (LM.GT.MXMMCT) THEN
                CALL QUIT('MXMMCT too small for LCLOSE')
              ENDIF 
            ENDIF
  209     CONTINUE
  208   CONTINUE

        DO 210 I=1,LM
          IQM = WORK(KIQM+I-1)
          IMM = WORK(KIMM+I-1)
          XT  = MMCORD(1,IMM) 
          YT  = MMCORD(2,IMM)
          ZT  = MMCORD(3,IMM)

          DITEMP = DUMMY
          DO 211 J=1,MMCENT ! find nearest MM center NOT on the KIMM list
            IF (J.NE.IMM) THEN
              LFIND = .FALSE.
              DO 212 K=1,LM
                IF (J.EQ.WORK(KIMM+K-1)) LFIND = .TRUE.
C               In case we have provided atomic numbers for the MM sites move the MM properties to only 
C               nuclear sites         
                IF (NELEME .NE. 0) THEN
                  IF (ELEME(J) .EQ. 0) LFIND = .TRUE.
                ENDIF
 212          CONTINUE
              IF (.NOT. LFIND) THEN
                DIMM2 = (MMCORD(1,J)-XT)**2 +
     *                  (MMCORD(2,J)-YT)**2 +
     *                  (MMCORD(3,J)-ZT)**2  
                DIMM  = SQRT(DIMM2)
                IF (DIMM .LT. DITEMP) THEN
                  DITEMP = DIMM
                  JSAVE = J
                ENDIF
              ENDIF
            ENDIF
  211     CONTINUE

C         now move MM properties from MM site IMM to site JSAVE. put lmmnul = true
C         for the affecte sites.
          IF (NMULT .GE. 0) THEN
            MUL0MM(JSAVE) = MUL0MM(JSAVE) + MUL0MM(IMM)
            MUL0MM(IMM)   = 0.0D0
            LMMNUL(JSAVE) = .TRUE.
            LMMNUL(IMM)   = .TRUE.
          ENDIF
          IF (NMULT .GE. 1) THEN
            MUL1MM(1,JSAVE) = MUL1MM(1,JSAVE) + MUL1MM(1,IMM)
            MUL1MM(2,JSAVE) = MUL1MM(2,JSAVE) + MUL1MM(2,IMM)
            MUL1MM(3,JSAVE) = MUL1MM(3,JSAVE) + MUL1MM(3,IMM)
            MUL1MM(1,IMM)   = 0.0D0
            MUL1MM(2,IMM)   = 0.0D0
            MUL1MM(3,IMM)   = 0.0D0
            LMMNUL(JSAVE)   = .TRUE.
            LMMNUL(IMM)     = .TRUE.
          ENDIF
          IF (NMULT .GE. 2) THEN
            MUL2MM(1,JSAVE) = MUL2MM(1,JSAVE) + MUL2MM(1,IMM)
            MUL2MM(2,JSAVE) = MUL2MM(2,JSAVE) + MUL2MM(2,IMM)
            MUL2MM(3,JSAVE) = MUL2MM(3,JSAVE) + MUL2MM(3,IMM)
            MUL2MM(4,JSAVE) = MUL2MM(4,JSAVE) + MUL2MM(4,IMM)
            MUL2MM(5,JSAVE) = MUL2MM(5,JSAVE) + MUL2MM(5,IMM)
            MUL2MM(6,JSAVE) = MUL2MM(6,JSAVE) + MUL2MM(6,IMM)
            MUL2MM(1,IMM)   = 0.0D0
            MUL2MM(2,IMM)   = 0.0D0
            MUL2MM(3,IMM)   = 0.0D0
            MUL2MM(4,IMM)   = 0.0D0
            MUL2MM(5,IMM)   = 0.0D0
            MUL2MM(6,IMM)   = 0.0D0
            LMMNUL(JSAVE)   = .TRUE.
            LMMNUL(IMM)     = .TRUE.
          ENDIF
          IF (IPOLTP .EQ. 1) THEN
            POLIMM(JSAVE) = POLIMM(JSAVE) + POLIMM(IMM)
            POLIMM(IMM)   = 0.0D0
            LMMNUL(JSAVE) = .TRUE.
            LMMNUL(IMM)   = .TRUE.
          ELSE IF (IPOLTP .EQ. 2) THEN
            POLMM(1,JSAVE) = POLMM(1,JSAVE) + POLMM(1,IMM)
            POLMM(2,JSAVE) = POLMM(2,JSAVE) + POLMM(2,IMM)
            POLMM(3,JSAVE) = POLMM(3,JSAVE) + POLMM(3,IMM)
            POLMM(4,JSAVE) = POLMM(4,JSAVE) + POLMM(4,IMM)
            POLMM(5,JSAVE) = POLMM(5,JSAVE) + POLMM(5,IMM)
            POLMM(6,JSAVE) = POLMM(6,JSAVE) + POLMM(6,IMM)
            POLMM(1,IMM)   = 0.0D0
            POLMM(2,IMM)   = 0.0D0
            POLMM(3,IMM)   = 0.0D0
            POLMM(4,IMM)   = 0.0D0
            POLMM(5,IMM)   = 0.0D0
            POLMM(6,IMM)   = 0.0D0
            LMMNUL(JSAVE)  = .TRUE.
            LMMNUL(IMM)    = .TRUE.
          ENDIF
  210   CONTINUE

      ENDIF

      IF (LCLOSE) THEN

        IOPTMM = 2
        CALL OUTMM(IOPTMM)

        IF (IQMMMC .EQ. 1) THEN
          WRITE(LUPRI,*)
          WRITE(LUPRI,*) 'WARNING: Removed ',TEMPQ,' charge from MM'
          WRITE(LUPRI,*)
          WRITE(LUPRI,*) 'WARNING: No charge renormalization'
        ENDIF

C       print now the closest QM-MM distance for MM sites with non-zero parameters

        DIMIN = DUMMY
        ISAVE = -1000
        JSAVE = -1000
        DO 213 I = 1,NUCIND
          DO 214 J = 1,MMCENT
            DIST2 = (CORD(1,I)-MMCORD(1,J))**2 +
     *              (CORD(2,I)-MMCORD(2,J))**2 +
     *              (CORD(3,I)-MMCORD(3,J))**2
            DIST = SQRT(DIST2)

            LHIT = .FALSE.
            IF (NMULT .GE. 0) THEN
              TEMP = MUL0MM(J)*MUL0MM(J)
              IF (TEMP .GT. THRMM) LHIT = .TRUE.
            ENDIF
            IF (NMULT .GE. 1) THEN
              TEMP = MUL1MM(1,J)*MUL1MM(1,J) + 
     *               MUL1MM(2,J)*MUL1MM(2,J) + 
     *               MUL1MM(3,J)*MUL1MM(3,J)
              IF (TEMP .GT. THRMM) LHIT = .TRUE.
            ENDIF
            IF (NMULT .GE. 2) THEN
                TEMP = MUL2MM(1,J)*MUL2MM(1,J) +
     *                 MUL2MM(2,J)*MUL2MM(2,J) +
     *                 MUL2MM(3,J)*MUL2MM(3,J) +
     *                 MUL2MM(4,J)*MUL2MM(4,J) +
     *                 MUL2MM(5,J)*MUL2MM(5,J) +
     *                 MUL2MM(6,J)*MUL2MM(6,J) 
                IF (TEMP .GT. THRMM) LHIT = .TRUE.
            ENDIF
            IF (IPOLTP .EQ. 1) THEN
              TEMP = POLIMM(J)*POLIMM(J)
              IF (TEMP .GT. THRMM) LHIT = .TRUE.
            ELSE IF (IPOLTP .EQ. 2) THEN
              TEMP = POLMM(1,J)*POLMM(1,J) +
     *               POLMM(2,J)*POLMM(2,J) +
     *               POLMM(3,J)*POLMM(3,J) +
     *               POLMM(4,J)*POLMM(4,J) +
     *               POLMM(5,J)*POLMM(5,J) +
     *               POLMM(6,J)*POLMM(6,J) 
              IF (TEMP .GT. THRMM) LHIT = .TRUE. 
            ENDIF

            IF (LHIT) THEN
              IF (DIST .LT. DIMIN) THEN
                DIMIN = DIST
                ISAVE = I
                JSAVE = J
              ENDIF
            ENDIF 

  214     CONTINUE
  213   CONTINUE

      WRITE(LUPRI,*)
      WRITE(LUPRI,*) 'Closest QM-MM distance is now ', DIMIN, ' au',
     *               ' for QM-MM pair ', ISAVE, JSAVE

      ENDIF ! if lclose

c      DO 215 I = 1,MMCENT
c        DO 216 J = I,MMCENT
c          IF (I.NE.J) THEN
c            DIST2 = (MMCORD(1,I)-MMCORD(1,J))**2 +
c     *              (MMCORD(2,I)-MMCORD(2,J))**2 +
c     *              (MMCORD(3,I)-MMCORD(3,J))**2
c            DIST = SQRT(DIST2) 
c            IF ( (DIST .LE. 3.40) .AND. (DIST .GT. 1.88) ) THEN
c              WRITE(LUPRI,*) 'I,J,DIST ',I,J,DIST
c            ENDIF
c          ENDIF
c  216   CONTINUE
c  215 CONTINUE



4000    FORMAT(I4,3(2X,F12.6))
4001    FORMAT(I4,2X,F12.6)
4002    FORMAT(I4,6(2X,F12.6))

      CALL QEXIT('READMM')
      END
C*************************************************************************************
C  /* Deck outmm */
      SUBROUTINE OUTMM(IOPT)
C------------------------------------------------------------------------
C     IF IOPT = 1 all information is printed
C     otherwise only the sites changed are printed
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "mxcent.h"
#include "abainf.h"
#include "inftap.h"
#include "nuclei.h"
#include "qmmm.h"
#include "codata.h"
C
      INTEGER IOPT

      CALL QENTER('OUTMM')

C     Output information

      IF (IOPT .EQ. 1) THEN
        WRITE(LUPRI,*) 
        WRITE(LUPRI,*) ' ---------------- '
        WRITE(LUPRI,*) ' QMMM information '
        WRITE(LUPRI,*) ' ---------------- '
      ELSE
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) 'WARNING: The following sites have been changed'
      ENDIF

      WRITE(LUPRI,*)
      WRITE(LUPRI,*) ' MM coordinates in au '
      WRITE(LUPRI,*) ' -------------------- '

      IF (NELEME .EQ. 0) THEN
        DO 199 J = 1,MMCENT
          IF (IOPT .EQ. 1) THEN
            WRITE(LUPRI,4000) J,MMCORD(1,J),MMCORD(2,J),MMCORD(3,J)
          ELSE IF (LMMNUL(J)) THEN
              WRITE(LUPRI,4000) J,MMCORD(1,J),MMCORD(2,J),MMCORD(3,J)
          ENDIF
  199   CONTINUE
      ELSE
        DO 200 J = 1,MMCENT
          IF (IOPT .EQ. 1) THEN
            WRITE(LUPRI,3999) J,ELEME(J),MMCORD(1,J),MMCORD(2,J),
     *                        MMCORD(3,J)
          ELSE IF (LMMNUL(J)) THEN
              WRITE(LUPRI,3999) J,ELEME(J),MMCORD(1,J),MMCORD(2,J),
     *                          MMCORD(3,J)
          ENDIF
  200   CONTINUE
      ENDIF

      IF (NMULT .GE. 0) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' MM charges '
        WRITE(LUPRI,*) ' ---------- '
        DO 201 J = 1,MMCENT
          IF (IOPT .EQ. 1) THEN
            WRITE(LUPRI,4001) J, MUL0MM(J)
          ELSE IF (LMMNUL(J)) THEN
              WRITE(LUPRI,4001) J, MUL0MM(J)
          ENDIF
  201   CONTINUE
      ENDIF

      IF (NMULT .GE. 1) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' MM dipoles (x,y,z) '
        WRITE(LUPRI,*) ' ------------------ '
        DO 202 J = 1,MMCENT
          IF (IOPT .EQ. 1) THEN
            WRITE(LUPRI,4000) J,MUL1MM(1,J),MUL1MM(2,J),MUL1MM(3,J)
          ELSE IF (LMMNUL(J)) THEN
              WRITE(LUPRI,4000) J,MUL1MM(1,J),MUL1MM(2,J),MUL1MM(3,J)
          ENDIF 
  202   CONTINUE
      ENDIF

      IF (NMULT .GE. 2) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' MM quadrupoles (xx,xy,xz,yy,yz,zz) '
        WRITE(LUPRI,*) ' ---------------------------------- '
        DO 203 J = 1,MMCENT
          IF (IOPT .EQ. 1) THEN
            WRITE(LUPRI,4002) J,MUL2MM(1,J),MUL2MM(2,J),MUL2MM(3,J),
     *                          MUL2MM(4,J),MUL2MM(5,J),MUL2MM(6,J)
          ELSE IF (LMMNUL(J)) THEN
              WRITE(LUPRI,4002) J,MUL2MM(1,J),MUL2MM(2,J),MUL2MM(3,J),
     *                            MUL2MM(4,J),MUL2MM(5,J),MUL2MM(6,J)
          ENDIF
  203   CONTINUE  
      ENDIF

      IF (IPOLTP .EQ. 1) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' MM isotropic polarizabilities in au '
        WRITE(LUPRI,*) ' ----------------------------------- '
        DO 204 J = 1,MMCENT
          IF (IOPT .EQ. 1) THEN
            WRITE(LUPRI,4001) J, POLIMM(J)
          ELSE IF (LMMNUL(J)) THEN
              WRITE(LUPRI,4001) J, POLIMM(J)
          ENDIF 
  204   CONTINUE
      ENDIF

      IF (IPOLTP .EQ. 2) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' MM polarizabilities in au (xx,xy,xz,yy,yz,zz) '
        WRITE(LUPRI,*) ' --------------------------------------------- '
        DO 205 J = 1,MMCENT
          IF (IOPT .EQ. 1) THEN
            WRITE(LUPRI,4002) J,POLMM(1,J),POLMM(2,J),POLMM(3,J),
     *                          POLMM(4,J),POLMM(5,J),POLMM(6,J)
          ELSE IF (LMMNUL(J)) THEN
              WRITE(LUPRI,4002) J,POLMM(1,J),POLMM(2,J),POLMM(3,J),
     *                            POLMM(4,J),POLMM(5,J),POLMM(6,J)
          ENDIF
  205   CONTINUE
      ENDIF

3999    FORMAT(I6,I4,3(2X,F12.6))
4000    FORMAT(I6,3(2X,F12.6))
4001    FORMAT(I6,2X,F12.6)
4002    FORMAT(I6,6(2X,F12.6))
4004    FORMAT(A80)

      CALL QEXIT('OUTMM')
      END
C**** end of herqm3.F *****************************************************************




