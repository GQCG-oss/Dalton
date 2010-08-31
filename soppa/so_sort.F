C
C  /* Deck so_sort */
      SUBROUTINE SO_SORT(POINT,NR,VEC,LVEC,PARRAY)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C
C     PURPOSE: Find the NR lowest values in the vector VEC 
C              and leave pointers to them in PARRAY.
C
#include <implicit.h>
#include <priunit.h>
C
      INTEGER   POINT(NR), PARRAY(LVEC)
      DIMENSION VEC(LVEC)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_SORT')
C
C------------------------------
C     Initialize pointer-array.
C------------------------------
C
      DO 100 I = 1,LVEC
         PARRAY(I) = I
  100 CONTINUE
C
C-------------------------------------------
C     Find pointers to the lowest NR values.
C-------------------------------------------
C
      DO 200 J = 1,NR
C
         DO 300 I = LVEC,1+J,-1
C
            IF ( VEC(PARRAY(I)) .LT. VEC(PARRAY(I-1)) ) THEN
C
               PI          = PARRAY(I-1)
               PARRAY(I-1) = PARRAY(I)
               PARRAY(I)   = PI
C
            END IF
C
  300    CONTINUE
C
  200 CONTINUE
C
C-----------------------------------------------------
C     Copy pointers for the lowest NR values to POINT.
C-----------------------------------------------------
C
      DO 400 I = 1,NR
         POINT(I) = PARRAY(I)
  400 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_SORT')
C
      RETURN
C
      END
