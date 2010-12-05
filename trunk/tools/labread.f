      PROGRAM READLABL
C
C Written by Hans Jorgen Aa. Jensen, Nov. 1983
C CRAY version 10. Oct. 1986
C Alliant version 30. Dec. 1987
C
      DOUBLE PRECISION TIME, SECOND
      CHARACTER*8 B(4), STARS
      PARAMETER (STARS = '********')
      PARAMETER (LUINP = 5, LUPRI = 6)
C
C      TIME = SECOND()
      WRITE (LUPRI, '(/A/)') ' MOLECULE labels found on file :'
C
      REWIND LUINP
      IREC = 0
      IERR = 0
    1 READ (LUINP,END=3,ERR=2) B
         IREC = IREC + 1
      IF (B(1) .NE. STARS) GO TO 1
         WRITE (LUPRI, '(5X,I5,3X,4(2X,A8))')  IREC, B
      GO TO 1
C
    2 CONTINUE
      IREC = IREC + 1
      IERR = IERR + 1
      WRITE (LUPRI, '(/A,I5/)') ' ERROR READING RECORD NO.',IREC
      REWIND LUINP
      DO 102 I = 1,IREC
         READ (LUINP) J
  102 CONTINUE
      IF (IERR .LE. 2) GO TO 1
  202 CONTINUE
         READ (LUINP,END=3) J
         IREC = IREC + 1
      GO TO 202
C
    3 CONTINUE
C      TIME = SECOND() - TIME
      WRITE (LUPRI,'(/I10,A)') IREC,
     *   ' records read before EOF on file.'
      WRITE (LUPRI,'(F10.2,A)') TIME,' CPU seconds used.'
      STOP ' '
C
      END
