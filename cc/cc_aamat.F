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
c /* deck CC_AAMAT */
*=====================================================================*
      SUBROUTINE CC_AAMAT(IATRAN,NATRAN,LISTO,LISTR,IOPTRES,FILAAM,
     &                     WORK, LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: driver for a list of A{O} matrix transformations
*
*           IOPTRES = 3 :  each result vector is written to its own
*                          file by a call to CC_WRRSP, FILAAM is used
*                          as list type and IATRAN(3,*) as list index
*
*           IOPTRES = 4 :  each result vector is added to a vector on
*                          file by a call to CC_WARSP, FILAAM is used
*                          as list type and IATRAN(3,*) as list index
*
*    N.B.: this is a very first quick hack solution...
*          should finaly be incorporated in the CCCR_AA routine...
*
*    Written by Christof Haettig, maj 1996.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#  include <priunit_dec.h>
#  include <ccorb_dec.h>
#  include <ccsdsym_dec.h>
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccroper.h>
#include <cclists.h>

* local parameters:
      CHARACTER*(19) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_AAMAT> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
    
      INTEGER LWORK, NATRAN

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL ZERO
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ZERO
#endif
      PARAMETER (ZERO = 0.0d0)

      CHARACTER LISTO*(*), LISTR*(*), FILAAM*(*)
      CHARACTER*(10) MODEL
      CHARACTER*(8)  LABEL
      LOGICAL LTWOEL
      INTEGER IATRAN(MXDIM_AATRAN,NATRAN)
      INTEGER IOPER, ITAMP, IFILE, IVEC, ITRAN, IOPT, IOPTRES
      INTEGER ISYOPE, ISYTAM, ISYRES, IRELAX
      INTEGER KEND1, LEND1, KRES0, KRES1, KRES2, LEN


* external functions:
      INTEGER ILSTSYM
#if defined (SYS_CRAY)
      REAL DDOT
#else
      DOUBLE PRECISION DDOT
#endif

*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
* check list combinations:
      IF ( (LISTO(1:2).NE.'o1' .AND. LISTO(1:2).NE.'o2')
     &    .OR. ( LISTR(1:1).NE.'R'.AND. LISTR(1:2).NE.'ER')
     &    .OR. (FILAAM(1:1).NE.'O'.AND.FILAAM(1:2).NE.'EO'
     &                            .AND.FILAAM(1:2).NE.'CO') ) THEN
        WRITE (LUPRI,*) 
     &        'A{O} matrix driver called with a strange list set:'
        WRITE (LUPRI,*) 'LISTO:',LISTO
        WRITE (LUPRI,*) 'LISTR:',LISTR
        WRITE (LUPRI,*) 'FILAAM :',FILAAM
        CALL QUIT('A{O} mat. driver called with a strange '//
     &             'set of list types.')
      END IF

      IF ( LISTO(1:2).EQ.'O2' .AND. NATRAN.NE.0) THEN
         WRITE (LUPRI,*) 
     &        'LISTO="o2" case not yet implemented in CC_AAMAT.'
         CALL QUIT('LISTO="o2" case not yet implemented in CC_AAMAT.')
      END IF

*---------------------------------------------------------------------*
* start loop over all requested A{O} matrix transformations
*---------------------------------------------------------------------*
      DO ITRAN = 1, NATRAN

        IOPER  = IATRAN(1,ITRAN)
        ITAMP  = IATRAN(2,ITRAN)
        IFILE  = IATRAN(3,ITRAN)
        IRELAX = IATRAN(4,ITRAN)

        LABEL  = LBLOPR(IOPER)
        LTWOEL = LPDBSOP(IOPER)

        ISYOPE = ILSTSYM(LISTO,IOPER)
        ISYTAM = ILSTSYM(LISTR,ITAMP)
        ISYRES = MULD2H(ISYOPE,ISYTAM)

        IF (ILSTSYM(FILAAM,IFILE).NE.ISYRES) THEN
          CALL QUIT('Symmetry mismatch in CC_AAMAT.')
        END IF

*---------------------------------------------------------------------*
* calculate A{O} matrix times a response amplitude vector:
*---------------------------------------------------------------------*
        KRES1 = 1
        KRES2 = KRES1 + NT1AM(ISYRES)
        KEND1 = KRES2 + NT2AM(ISYRES)
        LEND1 = LWORK - KEND1

        IF (LEND1.LT.0) CALL QUIT('Insufficient memory in CC_AAMAT.')

        IF ( (IRELAX.GT.0) .OR. LTWOEL ) THEN

           CALL CC_FDAAMAT(LISTR,ITAMP,WORK(KRES1),WORK(KRES2),
     &                     LABEL,IRELAX,WORK(KEND1),LEND1)

        ELSE

           CALL CCCR_AA(LABEL,ISYOPE,LISTR,ITAMP,WORK,LWORK)

           CALL CCLR_DIASCL(WORK(KRES2),2.0d0,ISYRES)

        END IF

*---------------------------------------------------------------------*
* write/add to vector on file:
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
           CALL QUIT('Unkown coupled cluster model in CC_AAMAT.')
        END IF

        KRES0 = -999999

        IF      (IOPTRES.EQ.3) THEN
          CALL CC_WRRSP(FILAAM, IFILE, ISYRES, IOPT, MODEL,
     &                  WORK(KRES0),WORK(KRES1),WORK(KRES2),
     &                  WORK(KEND1),LEND1)
        ELSE IF (IOPTRES.EQ.4) THEN
          CALL CC_WARSP(FILAAM, IFILE, ISYRES, IOPT, MODEL,
     &                  WORK(KRES0),WORK(KRES1),WORK(KRES2),
     &                  WORK(KEND1),LEND1)
        ELSE
          WRITE (LUPRI,*) 'Illegal option IOPT in CC_AAMAT.'
          CALL QUIT('Illegal option IOPT in CC_AAMAT.')
        END IF
        
        IF (LOCDBG) THEN
         WRITE (LUPRI,*) MSGDBG, 'wrote ',FILAAM,':',IFILE,' to disk.'
         LEN = NT1AM(ISYRES)
         IF (.NOT.CCS) LEN = LEN + NT2AM(ISYRES)
         WRITE (LUPRI,*) MSGDBG, 'NORM^2 = ', 
     &    DDOT(LEN,WORK(KRES1),1,WORK(KRES1),1)
        END IF

*---------------------------------------------------------------------*
* end of loop over A{O} transformations:
*---------------------------------------------------------------------*
      END DO

      RETURN
      END
*---------------------------------------------------------------------*
*               END OF SUBROUTINE CC_AAMAT                            *
*---------------------------------------------------------------------*
