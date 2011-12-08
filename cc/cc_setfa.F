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
c /* deck cc_setfa */
*=====================================================================*
      SUBROUTINE CC_SETFA(IFTRAN,IFDOTS,MXTRAN,MXVEC,IZETAV,
     &                    IOPER,IKAPPA,ITAMPA,ITAMPB,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: set up list of F matrix transformations 
*
*             IFTRAN - list of F matrix transformations
*             IFDOTS - list of vectors it should be dottet on
*        
*             MXTRAN - maximum list dimension
*             MXVEC  - maximum second dimension for IFDOTS
*      
*             IZETAV - index of lagrangian multiplier vector
*             IOPER  - index of property operator 
*             IKAPPA - index of the relaxation vector
*             ITAMPA - index of amplitude vector A
*             ITAMPB - index of amplitude vector B
*
*             ITRAN - index in IFTRAN list
*             IVEC  - second index in IFDOTS list
*
*    Written by Christof Haettig, november 1996.
*    IKAPPA entry added in june 1999
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>

      INTEGER MXVEC, MXTRAN
      INTEGER IFTRAN(5,MXTRAN)
      INTEGER IFDOTS(MXVEC,MXTRAN)

      LOGICAL LFNDA, LFNDB
      INTEGER IZETAV, IOPER, IKAPPA, ITAMPA, ITAMPB
      INTEGER ITRAN, IVEC
      INTEGER ITAMP, I, IDX

* statement  functions:
      LOGICAL LFATST, LFAEND
      INTEGER IL, IA, IO, IK
      LFATST(ITRAN,IL,IO,IK,IA) = 
     &        IFTRAN(1,ITRAN).EQ.IL .AND. IFTRAN(2,ITRAN).EQ.IO 
     &  .AND. IFTRAN(3,ITRAN).EQ.IA .AND. IFTRAN(5,ITRAN).EQ.IK
      LFAEND(ITRAN) = ITRAN.GT.MXTRAN .OR.
     &      (IFTRAN(1,ITRAN)+IFTRAN(2,ITRAN)+IFTRAN(3,ITRAN)).LE.0 


*---------------------------------------------------------------------*
* set up list of F{A} matrix transformations
*---------------------------------------------------------------------*
      ITRAN = 1
      LFNDA  = LFATST(ITRAN,IZETAV,IOPER,IKAPPA,ITAMPB)
      LFNDB  = LFATST(ITRAN,IZETAV,IOPER,IKAPPA,ITAMPA)

      DO WHILE ( .NOT. (LFNDA.OR.LFNDB.OR.LFAEND(ITRAN)))
       ITRAN = ITRAN + 1
       LFNDA  = LFATST(ITRAN,IZETAV,IOPER,IKAPPA,ITAMPB)
       LFNDB  = LFATST(ITRAN,IZETAV,IOPER,IKAPPA,ITAMPA)
      END DO

      IF (.NOT.(LFNDA.OR.LFNDB)) THEN
        IFTRAN(1,ITRAN) = IZETAV
        IFTRAN(2,ITRAN) = IOPER
        IFTRAN(3,ITRAN) = ITAMPA
        IFTRAN(4,ITRAN) = 0
        IFTRAN(5,ITRAN) = IKAPPA
        ITAMP = ITAMPB
      ELSE 
        IF (LFNDA) ITAMP = ITAMPA
        IF (LFNDB) ITAMP = ITAMPB
      END IF

      IVEC = 1
      DO WHILE (IFDOTS(IVEC,ITRAN).NE.ITAMP .AND.
     &            IFDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
        IVEC = IVEC + 1
      END DO

      IFDOTS(IVEC,ITRAN) = ITAMP

*---------------------------------------------------------------------*
      IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
        WRITE (LUPRI,*) 'Overflow error in CC_SETFA:'
        WRITE (LUPRI,*) 'IVEC, MXVEC  :',IVEC, MXVEC
        WRITE (LUPRI,*) 'ITRAN, MXTRAN:',ITRAN, MXTRAN
        WRITE (LUPRI,*) 'IOPER,IKAPPA :',IOPER,IKAPPA
        WRITE (LUPRI,*) 'IZETAV,ITAMPA,ITAMPB:',IZETAV,ITAMPA,ITAMPB
        IDX = 1
        DO WHILE ( .NOT. LFAEND(IDX) )
          WRITE(LUPRI,'(A,5I5,5X,(12I5,20X))') 'CC_SETFA>',
     &       (IFTRAN(I,IDX),I=1,5),(IFDOTS(I,IDX),I=1,MXVEC)
          IDX = IDX + 1
        END DO
        CALL FLSHFO(LUPRI)
        CALL QUIT('Overflow error in CC_SETFA')
      END IF
      
      RETURN
      END 

*---------------------------------------------------------------------*
*                END OF SUBROUTINE CC_SETFA                           *
*---------------------------------------------------------------------*
