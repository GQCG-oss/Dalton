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
c /* deck CC_SETXE */
*=====================================================================*
       SUBROUTINE CC_SETXE(TYPE,IXETRAN,IDOTS,MXTRAN,MXVEC,
     &                     IZETA,IOPER,IRLX1,IRLX2,IRLX3,IRLX4,
     &                     IDOTVEC,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: set up list of Xi and Eta vectors 
*
*             IXETRAN - list for CC_XIETA routine
*             IDOTS   - list of vectors Xi/Eta should be dotted on
*        
*             MXTRAN  - maximum IXETRAN list dimension
*             MXVEC   - second maximum dimension maximum for IDOTS
*      
*             IZETA   - index of left vector for ETA (ignored for Xi)
*             IOPER   - index of operator
*             IRLX1   - index for 1. kappa vector (0 for unrelaxed)
*             IRLX3   - index for 2. kappa vector (0 for unrelaxed)
*             IRLX3   - index for 3. kappa vector (0 for unrelaxed)
*             IRLX4   - index for 4. kappa vector (0 for unrelaxed)
*             IDOTVEC - index of vector Xi/Eta should be dotted on
*
*             ITRAN   - index in IXETRAN list
*             IVEC    - second index in IDOTS list
*
*    Written by Christof Haettig, june 1999.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <cclists.h>
      INTEGER MXVEC, MXTRAN
      INTEGER IXETRAN(MXDIM_XEVEC,MXTRAN)
      INTEGER IDOTS(MXVEC,MXTRAN)

      LOGICAL LFND
      CHARACTER*3 TYPE
      INTEGER IOPER, IZETA, IRLX1, IRLX2, IRLX3, IRLX4, IDOTVEC
      INTEGER ITRAN, IVEC
      INTEGER I, IDX

* statement  functions:
      LOGICAL LXITST, LETATST, LXEEND
      INTEGER IL, IO, IR1, IR2, IR3, IR4

      LXITST(ITRAN,IO,IR1,IR2,IR3,IR4) = IXETRAN(1,ITRAN).EQ.IO 
     & .AND. (
     &    (  ( (IXETRAN(5,ITRAN).EQ.IR1 .AND. IXETRAN(6,ITRAN).EQ.IR2)
     &     .OR.(IXETRAN(5,ITRAN).EQ.IR2 .AND. IXETRAN(6,ITRAN).EQ.IR1))
     & .AND. ( (IXETRAN(7,ITRAN).EQ.IR3 .AND. IXETRAN(8,ITRAN).EQ.IR4)
     &     .OR.(IXETRAN(7,ITRAN).EQ.IR4 .AND. IXETRAN(8,ITRAN).EQ.IR3))
     &    ) .OR.
     &    (  ( (IXETRAN(5,ITRAN).EQ.IR1 .AND. IXETRAN(6,ITRAN).EQ.IR3)
     &     .OR.(IXETRAN(5,ITRAN).EQ.IR3 .AND. IXETRAN(6,ITRAN).EQ.IR1))
     & .AND. ( (IXETRAN(7,ITRAN).EQ.IR2 .AND. IXETRAN(8,ITRAN).EQ.IR4)
     &     .OR.(IXETRAN(7,ITRAN).EQ.IR4 .AND. IXETRAN(8,ITRAN).EQ.IR2))
     &    ) .OR.
     &    (  ( (IXETRAN(5,ITRAN).EQ.IR1 .AND. IXETRAN(6,ITRAN).EQ.IR4)
     &     .OR.(IXETRAN(5,ITRAN).EQ.IR4 .AND. IXETRAN(6,ITRAN).EQ.IR1))
     & .AND. ( (IXETRAN(7,ITRAN).EQ.IR3 .AND. IXETRAN(8,ITRAN).EQ.IR2)
     &     .OR.(IXETRAN(7,ITRAN).EQ.IR2 .AND. IXETRAN(8,ITRAN).EQ.IR3))
     &    )  )

      LETATST(ITRAN,IL,IO,IR1,IR2,IR3,IR4) = 
     &   IXETRAN(2,ITRAN).EQ.IL .AND.  LXITST(ITRAN,IO,IR1,IR2,IR3,IR4)

      LXEEND(ITRAN) = ITRAN.GT.MXTRAN .OR. IXETRAN(1,ITRAN).LE.0

*---------------------------------------------------------------------*
* maintain list of Xi{A} and ETA{A} vectors:
*---------------------------------------------------------------------*
      IF (TYPE(1:3).EQ.'Xi ') THEN

        ITRAN = 1
        LFND  = LXITST(ITRAN,IOPER,IRLX1,IRLX2,IRLX3,IRLX4)

        DO WHILE ( .NOT.(LFND.OR.LXEEND(ITRAN))) 
         ITRAN = ITRAN + 1
         LFND  = LXITST(ITRAN,IOPER,IRLX1,IRLX2,IRLX3,IRLX4)
        END DO

        IF (.NOT.LFND) THEN
          IXETRAN(1,ITRAN) = IOPER
          IXETRAN(3,ITRAN) = 0
          IXETRAN(5,ITRAN) = IRLX1
          IXETRAN(6,ITRAN) = IRLX2
          IXETRAN(7,ITRAN) = IRLX3
          IXETRAN(8,ITRAN) = IRLX4
        END IF

      ELSE IF (TYPE(1:3).EQ.'Eta') THEN
        ITRAN = 1
        LFND  = LETATST(ITRAN,IZETA,IOPER,IRLX1,IRLX2,IRLX3,IRLX4)

        DO WHILE ( .NOT.(LFND.OR.LXEEND(ITRAN))) 
         ITRAN = ITRAN + 1
         LFND  = LETATST(ITRAN,IZETA,IOPER,IRLX1,IRLX2,IRLX3,IRLX4)
        END DO

        IF (.NOT.LFND) THEN
          IXETRAN(1,ITRAN) = IOPER
          IXETRAN(2,ITRAN) = IZETA
          IXETRAN(5,ITRAN) = IRLX1
          IXETRAN(6,ITRAN) = IRLX2
          IXETRAN(7,ITRAN) = IRLX3
          IXETRAN(8,ITRAN) = IRLX4
        END IF

      END IF


      IVEC = 1
      DO WHILE (IDOTS(IVEC,ITRAN).NE.IDOTVEC .AND.
     &            IDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
        IVEC = IVEC + 1
      END DO

      IDOTS(IVEC,ITRAN) = IDOTVEC
      IF (TYPE(1:3).EQ.'Eta') IXETRAN(4,ITRAN) = 0
      IF (TYPE(1:3).EQ.'Xi ') IXETRAN(3,ITRAN) = 0

*---------------------------------------------------------------------*
      IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
        WRITE (LUPRI,*) 'TYPE   :',TYPE
        WRITE (LUPRI,*) 'IZETA  :',IZETA
        WRITE (LUPRI,*) 'IOPER  :',IOPER 
        WRITE (LUPRI,*) 'IRLX1-4:',IRLX1,IRLX2,IRLX3,IRLX4
        WRITE (LUPRI,*) 'IDOTVEC:',IDOTVEC
        WRITE (LUPRI,*) 'IVEC   :',IVEC
        WRITE (LUPRI,*) 'ITRAN  :',ITRAN
        IDX = 1
        DO WHILE( .NOT. LXEEND(IDX) )
          WRITE(LUPRI,'(A,8I5,5X,(12I5,20X))') 'CC_SETXE>',
     &       (IXETRAN(I,IDX),I=1,8),(IDOTS(I,IDX),I=1,MXVEC)
          IDX = IDX + 1
        END DO
        CALL FLSHFO(LUPRI)
        CALL QUIT('Overflow error in CC_SETXE.')
      END IF
      
      RETURN
      END 

*---------------------------------------------------------------------*
*                END OF SUBROUTINE CC_SETXE                           *
*---------------------------------------------------------------------*
