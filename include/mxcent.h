C
C     file: mxcent.h
C
C     MXCENT = max number of nuclei + point charges + ghost orbital centers
C
C     IF you change MXCENT you should do a "make depend"
C     and then rebuild the program using the command "make".
C
      INTEGER MXCENT, MXCOOR
      PARAMETER (MXCENT = 120, MXCOOR = 3*MXCENT)
