*=====================================================================*
       SUBROUTINE CC_SETK11(IKTRAN,IKDOTS,MXTRAN,MXVEC,
     &                      IZETVA,IZETVB,ITAMPC,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: set up list of K matrix transformations
*             assumes that bar{T^A} and bar{T^B} members of the same list
*
*             IKTRAN - list of K matrix transformations
*             IKDOTS - list of vectors it should be dottet on
*
*             MXTRAN - maximum list dimension
*             MXVEC  - maximum second dimension for IBDOTS
*
*             IZETVA - index of lagrangian multiplier vector A
*             IZETVB - index of lagrangian multiplier vector B
*             ITAMPC - index of amplitude vector C
*
*             ITRAN - index in IKTRAN list
*             IVEC  - second index in IKDOTS list
*
*    CCMM JK+OC, modified version of CC_SETB11
*
*=====================================================================*
      IMPLICIT NONE

#include <priunit.h>
      INTEGER MXVEC, MXTRAN
      INTEGER IKTRAN(3,MXTRAN)
      INTEGER IKDOTS(MXVEC,MXTRAN)

      LOGICAL LFND
      INTEGER IZETVA, IZETVB, ITAMPC
      INTEGER ITRAN, IVEC
      INTEGER IZETV, I, IDX

* statement  functions:
      LOGICAL LBTST, LBEND
      INTEGER IB, IA
      LBTST(ITRAN,IA,IB) =
     &      ( IKTRAN(1,ITRAN).EQ.IA .AND. IKTRAN(2,ITRAN).EQ.IB ) .OR.
     &      ( IKTRAN(1,ITRAN).EQ.IB .AND. IKTRAN(2,ITRAN).EQ.IA )
      LBEND(ITRAN) = ITRAN.GT.MXTRAN .OR.
     &   (IKTRAN(1,ITRAN)+IKTRAN(2,ITRAN)).LE.0

*---------------------------------------------------------------------*
* set up list of K matrix transformations
*---------------------------------------------------------------------*
      ITRAN = 1
      LFND = LBTST(ITRAN,IZETVA,IZETVB)

      DO WHILE ( .NOT. (LFND.OR.LBEND(ITRAN)) )
       ITRAN = ITRAN + 1
       LFND = LBTST(ITRAN,IZETVA,IZETVB)
      END DO

      IF (.NOT.LFND) THEN
        IKTRAN(1,ITRAN) = IZETVA 
        IKTRAN(2,ITRAN) = IZETVB 
        IKTRAN(3,ITRAN) = 0
        IZETV = ITAMPC 
      ELSE
        IF (LFND) IZETV = ITAMPC 
      END IF

      IVEC = 1
      DO WHILE (IKDOTS(IVEC,ITRAN).NE.IZETV .AND.
     &           IKDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
        IVEC = IVEC + 1
      END DO

      IKDOTS(IVEC,ITRAN) = IZETV 

*---------------------------------------------------------------------*
      IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
        WRITE (LUPRI,*) 'IVEC :',IVEC
        WRITE (LUPRI,*) 'ITRAN:',ITRAN
        WRITE (LUPRI,*) 'IZETVA,IZETVB:',IZETVA,IZETVB
        IDX = 1
        DO WHILE ( .NOT. LBEND(IDX) )
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') 'CC_SETK11>',
     &       (IKTRAN(I,IDX),I=1,3),(IKDOTS(I,IDX),I=1,MXVEC)
          IDX = IDX + 1
        END DO
        CALL FLSHFO(LUPRI)
        CALL QUIT('Overflow error in CC_SETK11.')
      END IF

      RETURN
      END


