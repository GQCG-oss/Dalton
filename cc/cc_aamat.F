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
c /* deck CC_AAMAT */
*=====================================================================*
      SUBROUTINE CC_AAMAT(IATRAN,NATRAN,LISTO,LISTR,IOPTRES,FILAAM,
     &                    IADOTS,ACONS,MXVEC,WORK,LWORK)
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
*           IOPTRES = 5 : the result vectors are dotted on an array 
*                         of vectors, the type of the arrays is given 
*                         by FILAAM and the indices are taken from 
*                         IADOTS. the results are returned in ACONS
*
*    N.B.: this is a very first quick hack solution...
*          should finaly be incorporated in the CCCR_AA routine...
*
*    Written by Christof Haettig, maj 1996.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccroper.h>
#include <cclists.h>
#include <ccnoddy.h>
#include <dummy.h>

      CHARACTER*1 CDUMMY
      PARAMETER (CDUMMY = ' ')

* local parameters:
      CHARACTER*(19) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_AAMAT> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
    
      INTEGER LWORK, NATRAN, MXVEC

#if defined (SYS_CRAY)
      REAL WORK(LWORK), ACONS(MXVEC,NATRAN)
      REAL ZERO, TWO, FREQ
#else
      DOUBLE PRECISION WORK(LWORK), ACONS(MXVEC,NATRAN)
      DOUBLE PRECISION ZERO, TWO, FREQ
#endif
      PARAMETER (ZERO = 0.0d0, TWO = 2.0D0)

      CHARACTER LISTO*(*), LISTR*(*), FILAAM*(*)
      CHARACTER*(10) MODEL
      CHARACTER*(8)  LABEL
      LOGICAL LTWOEL
      INTEGER IATRAN(MXDIM_AATRAN,NATRAN), IADOTS(MXVEC,NATRAN)
      INTEGER IOPER, IDLSTR, IFILE, IVEC, ITRAN, IOPT, IOPTRES
      INTEGER ISYOPE, ISYTAM, ISYRES, IRELAX, IOPTE
      INTEGER KEND1, LWRK1, KRES1, KRES2, LEN, KRES1EFF, KRES2EFF


* external functions:
      INTEGER ILSTSYM
#if defined (SYS_CRAY)
      REAL DDOT, FREQLST
#else
      DOUBLE PRECISION DDOT, FREQLST
#endif

*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'Entered  cc_aamat...'
        CALL FLSHFO(LUPRI)
      END IF

      IF ( .NOT. (CCS .OR. CC2 .OR. CCSD .OR. CC3) ) THEN
        WRITE (LUPRI,'(/1x,a)') 'CC_AAMAT Called for a Coupled Cluster '
     &          //'method not implemented in CC_AAMAT...'
        CALL QUIT('Unkown CC method in CC_AAMAT.')
      END IF

* check list combinations:
      IF ( (LISTO(1:2).NE.'o1' .AND. LISTO(1:2).NE.'o2')
     &    .OR. ( LISTR(1:1).NE.'R'.AND. LISTR(1:2).NE.'ER')) THEN
        WRITE (LUPRI,*) 
     &        'A{O} matrix driver called with a strange list set:'
        WRITE (LUPRI,*) 'LISTO:',LISTO
        WRITE (LUPRI,*) 'LISTR:',LISTR
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
        IDLSTR = IATRAN(2,ITRAN)
        IFILE  = IATRAN(3,ITRAN)
        IRELAX = IATRAN(4,ITRAN)

        LABEL  = LBLOPR(IOPER)
        LTWOEL = LPDBSOP(IOPER)

        ISYOPE = ILSTSYM(LISTO,IOPER)
        ISYTAM = ILSTSYM(LISTR,IDLSTR)
        ISYRES = MULD2H(ISYOPE,ISYTAM)

        IF ((IOPTRES.EQ.3.OR.IOPTRES.EQ.4) .AND.
     &       ILSTSYM(FILAAM,IFILE).NE.ISYRES) THEN
          CALL QUIT('Symmetry mismatch in CC_AAMAT.')
        END IF

*---------------------------------------------------------------------*
* calculate A{O} matrix times a response amplitude vector:
*---------------------------------------------------------------------*

        ! allocate memory for results vector
        KRES1 = 1
        KRES2 = KRES1 + NT1AM(ISYRES)
        KEND1 = KRES2 + NT2AM(ISYRES)
        LWRK1 = LWORK - KEND1
        IF (LWRK1.LT.0) CALL QUIT('Insufficient memory in CC_AAMAT.')

        IF ( (IRELAX.GT.0) .OR. LTWOEL ) THEN

           CALL CC_FDAAMAT(LISTR,IDLSTR,WORK(KRES1),WORK(KRES2),
     &                     LABEL,IRELAX,WORK(KEND1),LWRK1)

           IF (CCSDT) CALL QUIT('relaxed with triples not available!')

        ELSE

c          --------------------------------------------------
c          do the CCS/CC2/CCSD parts:
c          --------------------------------------------------
           CALL CCCR_AA(LABEL,ISYOPE,LISTR,IDLSTR,DUMMY,WORK,LWORK)
           CALL CCLR_DIASCL(WORK(KRES2),2.0d0,ISYRES)

c          --------------------------------------------------
c          do the triples parts:
c          --------------------------------------------------
           IF (CCSDT) THEN
             ! find frequency associates with result vector
             IF (IOPTRES.EQ.5) THEN
               FREQ  = 0.0D0
             ELSE
               FREQ  = FREQLST(FILAAM,IFILE)
             END IF
          
             ! allocate memory for effective results vector
             KRES1EFF = KEND1
             KRES2EFF = KRES1EFF + NT1AM(ISYRES)
             KEND1    = KRES2EFF + NT2AM(ISYRES)
             LWRK1    = LWORK - KEND1
             IF (LWRK1.LT.0) 
     &          CALL QUIT('Insufficient memory in CC_AAMAT. (2)')

             CALL DZERO(WORK(KRES1EFF),NT1AM(ISYRES))
             CALL DZERO(WORK(KRES2EFF),NT2AM(ISYRES))

             IF (NODDY_AAMAT) THEN
               CALL CCSDT_AAMAT_NODDY(IOPTRES,FREQ,LABEL,ISYOPE,
     &                                LISTR,IDLSTR,.TRUE.,
     &                                WORK(KRES1),WORK(KRES2),
     &                                WORK(KRES1EFF),WORK(KRES2EFF),
     &                                IADOTS,ACONS,FILAAM,ITRAN,
     &                                NATRAN,MXVEC,WORK(KEND1),LWRK1)

             ELSE

               WRITE(LUPRI,*)'AAMATSD should be called from '
               WRITE(LUPRI,*)'CC_BMAT module. Check NEW_RHS flag '
               WRITE(LUPRI,*)'in CCRHSVEC.'
               CALL QUIT('Real code for AAMATSD not called here...')

             END IF


           END IF

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
        ELSE IF (CC3) THEN
           MODEL = 'CC3       '
           IOPT  = 3
           IOPTE = 24
        ELSE
           CALL QUIT('Unkown coupled cluster model in CC_AAMAT.')
        END IF

        IF      (IOPTRES.EQ.3) THEN
          CALL CC_WRRSP(FILAAM, IFILE, ISYRES, IOPT, MODEL,
     &                  DUMMY,WORK(KRES1),WORK(KRES2),
     &                  WORK(KEND1),LWRK1)
          IF (CCSDT) THEN
            CALL CC_WRRSP(FILAAM,IFILE,ISYRES,IOPTE,MODEL,DUMMY,
     &                 WORK(KRES1EFF),WORK(KRES2EFF),WORK(KEND1),LWRK1)
          END IF
        ELSE IF (IOPTRES.EQ.4) THEN
          CALL CC_WARSP(FILAAM, IFILE, ISYRES, IOPT, MODEL,
     &                  DUMMY,WORK(KRES1),WORK(KRES2),
     &                  WORK(KEND1),LWRK1)
          IF (CCSDT) THEN
            CALL CC_WARSP(FILAAM,IFILE,ISYRES,IOPTE,MODEL,DUMMY,
     &                 WORK(KRES1EFF),WORK(KRES2EFF),WORK(KEND1),LWRK1)
          END IF
        ELSE IF (IOPTRES.EQ.5) THEN
           IF (LOCDBG) THEN
             IVEC = 1
             WRITE(LUPRI,*) 'ACONS TRIPLES CONTRIBUTION:'
             DO WHILE (IADOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
                WRITE (LUPRI,*)
     &                'ACONS:',IVEC,ITRAN,ACONS(IVEC,ITRAN),IOPT
                IVEC = IVEC + 1
             END DO
           END IF
 
           CALL CCDOTRSP(IADOTS,ACONS,IOPT,FILAAM,ITRAN,NATRAN,MXVEC,
     &                   WORK(KRES1),WORK(KRES2),ISYRES,
     &                   WORK(KEND1),LWRK1)
 
           IF (LOCDBG) THEN
             IVEC = 1
             DO WHILE (IADOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
                WRITE (LUPRI,*)
     &                'ACONS:',IVEC,ITRAN,ACONS(IVEC,ITRAN),IOPT
                IVEC = IVEC + 1
             END DO
           END IF
        ELSE
          WRITE (LUPRI,*) 'Illegal option IOPTRES in CC_AAMAT.'
          CALL QUIT('Illegal option IOPTRES in CC_AAMAT.')
        END IF
        
        IF (LOCDBG .AND. (IOPTRES.EQ.3 .OR. IOPTRES.EQ.4)) THEN
         WRITE (LUPRI,*) MSGDBG, 'wrote ',FILAAM,':',IFILE,' to disk.'
         LEN = NT1AM(ISYRES)
         IF (.NOT.CCS) LEN = LEN + NT2AM(ISYRES)
         WRITE (LUPRI,*) MSGDBG, 'NORM^2 for result vector = ', 
     &    DDOT(LEN,WORK(KRES1),1,WORK(KRES1),1)
         WRITE (LUPRI,*) 'Listing of the result vector:'
         CALL CC_PRP(WORK(KRES1),WORK(KRES2),ISYRES,1,1)
         IF (CCSDT) THEN
           WRITE (LUPRI,*) 'norm^2 of eff. result vector:',
     &        DDOT(LEN,WORK(KRES1EFF),1,WORK(KRES1EFF),1)
           WRITE (LUPRI,*) 'Listing of the eff. result vector:'
           CALL CC_PRP(WORK(KRES1EFF),WORK(KRES2EFF),ISYRES,1,1)
         END IF
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
