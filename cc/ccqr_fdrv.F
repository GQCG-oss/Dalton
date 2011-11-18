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
c /* deck ccqr_fdrv */
*=====================================================================*
      SUBROUTINE CCQR_FDRV(LISTL, LISTB, LISTC,
     &                     NFTRAN, MXVEC,
     &                     IFTRAN, IFDOTS, FCON,
     &                     WORK, LWORK, OPTION)
*---------------------------------------------------------------------*
*
*    Purpose: driver for a list of F matrix transformations
*
*    Written by Christof Haettig, November 1996.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccroper.h>
#include <dummy.h>

* local parameters:
      CHARACTER*(19) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCQR_FDRV> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER MXTRAN, MXDOT
      PARAMETER (MXTRAN = 240, MXDOT = 60)
    

      INTEGER LWORK, NFTRAN, MXVEC

#if defined (SYS_CRAY)
      REAL FCON(MXVEC,NFTRAN)
      REAL WORK(LWORK)
      REAL ZERO
      REAL FCON1, FCON2
      REAL DDOT
      REAL SECOND, FTIME, TTIME
#else
      DOUBLE PRECISION FCON(MXVEC,NFTRAN)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ZERO
      DOUBLE PRECISION FCON1, FCON2
      DOUBLE PRECISION DDOT
      DOUBLE PRECISION SECOND, FTIME, TTIME
#endif
      PARAMETER (ZERO = 0.0d0)

      CHARACTER LISTL*(*), LISTB*(*), LISTC*(*), OPTION*(*)
      CHARACTER*(3)  LSTBTRN
      CHARACTER*(10) MODEL
      LOGICAL L_DO_FTRAN
      INTEGER IFTRAN(3,NFTRAN)
      INTEGER IFDOTS(MXVEC,NFTRAN)
      INTEGER ISYRES, ISYMTB, ISYMTC, ISYCTR, IFILE
      INTEGER IVEC, ITRAN, ITAMPB, ITAMPC, IZETAV, IOPT
      INTEGER KGAMMA1, KGAMMA2, KT1AMPC, KT2AMPC, KEND1, LEND1
      INTEGER KCTR1, KCTR2, KTHETA1, KTHETA2, LEN, IDXBTRN

* external functions:
      INTEGER ILSTSYM
      INTEGER IR2TAMP

*---------------------------------------------------------------------*
* start:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG, 'LISTL:',LISTL(1:3)
        WRITE (LUPRI,*) MSGDBG, 'LISTB:',LISTB(1:3)
        WRITE (LUPRI,*) MSGDBG, 'LISTC:',LISTC(1:3)
        WRITE (LUPRI,*) MSGDBG, 'NFTRAN:',NFTRAN
      END IF

      TTIME = SECOND()

      IF (IPRINT.GT.0) THEN
      
        WRITE (LUPRI,'(//1X,A1,50("="),A1)') '+','+'

        IF (LISTL(1:2).EQ.'L0') THEN
          WRITE (LUPRI,'(1X,A52)')
     &      '|        F MATRIX TRANSFORMATION SECTION           |'
        ELSE
          WRITE (LUPRI,'(1X,A52)')
     &      '|    GENERALIZED F MATRIX TRANSFORMATION SECTION   |'
        END IF

        IF (OPTION(1:7).EQ.'FILENEW') THEN
          WRITE (LUPRI,'(1X,A52)')
     &    '|          (result is written to file)             |'
        ELSE IF (OPTION(1:4).EQ.'FILE') THEN
          WRITE (LUPRI,'(1X,A52)')
     &    '|     (result is added to a vector on file)        |'
        END IF
        
        IF (OPTION(1:4).EQ.'FILE') THEN
          WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'
          WRITE (LUPRI,'(1X,A52)')
     &    '| L vector | R vector | result vector |            |'
          WRITE (LUPRI,'(1X,2A3,A,A3,A,A3,A)') '|  ',LISTL(1:3),
     &    ' No. |  ',LISTB(1:3),' No. |    ',LISTC(1:3),
     &    ' No.    | time/secs  |'
        ELSE IF (OPTION(1:4).EQ.'DOTP') THEN
          WRITE (LUPRI,'(1X,A52)')
     &    '|    (result used to calculate dot products)       |'
          WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'
          WRITE (LUPRI,'(1X,A52)')
     &    '| L vector | R vector | # dot prod.   |            |'
          WRITE (LUPRI,'(1X,2A3,A,A3,A,A3,A)') '|  ',LISTL(1:3),
     &    ' No. |  ',LISTB(1:3),' No. |   with ',LISTC(1:3),
     &    '    | time/secs  |'
        END IF

        WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'

      END IF
*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      IF (MXVEC*NFTRAN.NE.0) CALL DZERO(FCON,MXVEC*NFTRAN)

*---------------------------------------------------------------------*
* start loop over all requested F matrix transformations
*---------------------------------------------------------------------*
      DO ITRAN = 1, NFTRAN

        IZETAV = IFTRAN(1,ITRAN)
        ITAMPB = IFTRAN(2,ITRAN)
        IFILE  = IFTRAN(3,ITRAN)

        ISYCTR = ILSTSYM(LISTL,IZETAV)
        ISYMTB = ILSTSYM(LISTB,ITAMPB)

        ISYRES = MULD2H(ISYMTB,ISYCTR)

        KGAMMA1 = 1
        KGAMMA2 = KGAMMA1 + NT1AM(ISYRES)
        KT1AMPC = KGAMMA2 + NT2AM(ISYRES)
        KT2AMPC = KT1AMPC + NT1AM(ISYRES)
        KEND1   = KT2AMPC + NT2AM(ISYRES)
        LEND1   = LWORK - KEND1


        IF (LEND1 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CCQR_FDRV.')
        END IF

        FTIME = SECOND() 

*---------------------------------------------------------------------*
* check, if F matrix transformation is already stored on file:
* (disabled, for the time beeing, because of incompatibilities...)
*---------------------------------------------------------------------*
C      L_DO_FTRAN = .TRUE.

C      IF (OPTION(1:4).NE.'FILE') THEN

C         IF (LISTL.EQ.'L0' .AND. LISTB.EQ.'R1') THEN
C           IOPT = 33
C           CALL CC_RDRSP('F1',ITAMPB,ISYMTB,IOPT,MODEL,
C    &                    WORK(KGAMMA1),WORK(KGAMMA2))
C           IF (IOPT.NE.33) L_DO_FTRAN = .FALSE. ! F x R1 was on file
C         END IF

C      END IF
*---------------------------------------------------------------------*
* calculate F matrix times a response amplitude vector:
*      F * t^B = <Lambda^A| [[H,T^B],tau] |CC>
*---------------------------------------------------------------------*
C      IF (L_DO_FTRAN) THEN
C        IF (OPTION(1:4).EQ.'FILE') THEN
          WRITE (LUPRI,*) 'call CC_FMATOLD'
          CALL CC_FMATOLD(LISTL,IZETAV,LISTB,ITAMPB,WORK,LWORK)
C        ELSE
C         CALL CC_FTRAN(LISTL,IZETAV,LISTB,ITAMPB,WORK,LWORK)
C        END IF
C      END IF
c      IF (L_DO_FTRAN) THEN
c       WRITE (LUPRI,*) 'call CC_FGL'
c       CALL CC_FGL(LISTL,IZETAV,LISTB,ITAMPB,'--',0,.TRUE.,WORK,LWORK)
c      END IF

      IF (OPTION(1:4).EQ.'FILE') THEN
*---------------------------------------------------------------------*
* save transformed vector on file: 
*     default is, to add the vector to what is already on the file
*     only if OPTION='FILENEW' we initialize the file with the vector
*---------------------------------------------------------------------*
        IF (CCS) THEN
           MODEL = 'CCS       '
           IOPT  = 1
        ELSE IF (CC2) THEN
           MODEL = 'CC2       '
           IOPT  = 3
        ELSE IF (CCSD) THEN
           MODEL = 'CCSD      '
           IOPT  = 3
        ELSE
           CALL QUIT('Unknown coupled cluster model in CC_FDRV.')
        END IF

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) MSGDBG, 'LISTC:',LISTC
          WRITE (LUPRI,*) MSGDBG, 'IFILE:',IFILE
        END IF

        IF (OPTION(1:7).EQ.'FILENEW') THEN
          CALL CC_WRRSP(LISTC, IFILE, ISYRES, IOPT, MODEL,DUMMY,
     &                  WORK(KGAMMA1),WORK(KGAMMA2),WORK(KEND1),LEND1)
        ELSE
          CALL CC_WARSP(LISTC, IFILE, ISYRES, IOPT, MODEL,DUMMY,
     &                  WORK(KGAMMA1),WORK(KGAMMA2),WORK(KEND1),LEND1)
        END IF
        
        FTIME = SECOND() - FTIME

        IF (IPRINT.GT.0)  THEN
          WRITE (LUPRI,'(1X,A,I5,A,I5,A,I6,A,F10.2,A)') '| ',IZETAV,
     &      '    | ',ITAMPB,'    | ',IFILE,'        |',FTIME,'  |'
          CALL FLSHFO(LUPRI)
        END IF

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) MSGDBG, 'wrote ',LISTC,':',IFILE,' to disk.'
          IF (CCS) THEN 
            WRITE (LUPRI,*) MSGDBG, 'NORM^2 = ',
     &       DDOT(NT1AM(ISYRES),WORK(KGAMMA1),1,WORK(KGAMMA1),1)
          ELSE
            WRITE (LUPRI,*) MSGDBG, 'NORM^2 = ',
     &       DDOT(NT1AM(ISYRES)+NT2AM(ISYRES),WORK(KGAMMA1),1,
     &            WORK(KGAMMA1),1)
          END IF
        END IF

      ELSE IF (OPTION(1:4).EQ.'DOTP') THEN
*---------------------------------------------------------------------*
* calculate required dot products:
*---------------------------------------------------------------------*
       IVEC = 1
       DO WHILE (IFDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)

         ITAMPC = IFDOTS(IVEC,ITRAN)
         ISYMTC = ILSTSYM(LISTC,ITAMPC)
      
         IF (ISYMTC.NE.ISYRES) THEN
           CALL QUIT('symmetry mismatch in CCQR_FDRV.')
         END IF

         IF ((LISTC(1:1).EQ.'R').OR.(LISTC(1:1).EQ.'C')) THEN
            IOPT = 3
            CALL CC_RDRSP(LISTC,ITAMPC,ISYMTC,IOPT,MODEL,
     &                    WORK(KT1AMPC),WORK(KT2AMPC)   )
         ELSE IF (LISTC(1:1) .EQ. 'O') THEN
            CALL CC_XKSI(WORK(KT1AMPC),LBLOPR(ITAMPC),ISYMTC,0,
     *                   DUMMY,WORK(KEND1),LEND1)
         ELSE
            CALL QUIT(' Unknown list in CCQR_FDRV')
         ENDIF


         FCON1  = DDOT(NT1AM(ISYRES),WORK(KT1AMPC),1,WORK(KGAMMA1),1)
         IF (CC2 .OR. CCSD) THEN
          FCON2 = DDOT(NT2AM(ISYRES),WORK(KT2AMPC),1,WORK(KGAMMA2),1)
         ELSE IF (CCS) THEN
          FCON2 = ZERO
         ELSE
          CALL QUIT('ERROR> unknown CC model in CCQR_FDRV.')
         END IF

         FCON(IVEC,ITRAN) = FCON1 + FCON2

         IVEC = IVEC + 1
       END DO

       FTIME = SECOND() - FTIME

       IF (IPRINT.GT.0)  THEN
         WRITE (LUPRI,'(1X,A,I5,A,I5,A,I6,A,F10.2,A)') '| ',IZETAV,
     &     '    | ',ITAMPB,'    | ',IVEC-1,'        |',FTIME,'  |'
       END IF

*---------------------------------------------------------------------*
      ELSE
        CALL QUIT('Illegal OPTION in CCQR_FDRV')
      END IF

      END DO

*---------------------------------------------------------------------*
* print the results: 
*---------------------------------------------------------------------*

      TTIME = SECOND() - TTIME

      IF (IPRINT.GT.0) THEN
        WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'
        WRITE (LUPRI,'(1X,A,I4,A,F10.2,A)') 
     &   '| total time for',NFTRAN,' F transforms.:',TTIME,' secs.|'
        WRITE (LUPRI,'(1X,A1,50("="),A1,//)') '+','+'
      END IF

* F matrix transformations:
      IF (LOCDBG) THEN
        IF (MXVEC.NE.0) THEN
         WRITE (LUPRI,*) 'Results for F matrix transformations:'
         DO ITRAN = 1, NFTRAN
          WRITE (LUPRI,*)
          IVEC = 1
          DO WHILE (IFDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
           WRITE(LUPRI,'(A,3I5,2X,E18.12)') MSGDBG, 
     &      (IFTRAN(I,ITRAN),I=1,2),IFDOTS(IVEC,ITRAN),FCON(IVEC,ITRAN)
           IVEC = IVEC + 1
          END DO
         END DO
        END IF
        WRITE (LUPRI,*)
      END IF

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCQR_FDRV                            *
*---------------------------------------------------------------------*
