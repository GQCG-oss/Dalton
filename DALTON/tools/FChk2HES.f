C
C        ***** FChk2HES.f *****************
C
C     This program reads the Formatted Checkpointfile of Gaussian  
C     and creates a DALTON.HES file with Hessian and Atomic coordinates
C     as well as a dummy MOLECULE.INP file.
C     The checkpointfile "Test.FChk" is created with the keyword 
C     "FormCheck=ForceCart", and the keyword "NoSymm" has also to be specified.
C     Commandline: FChk2HES [-n] [filename] 
C                      -n does not read and write geometry;
C                      filename is optional, default is "Test.FChk"
C     *************************************************************************
C     G. Hangartner, 27.11.1996
C     *************************************************************************
C     REMARK:
C     Gaussian is a bit special:
C     1) If symmetry is used: Hessian will be in standard(Gaussian) orientation,
C        but geometry will be in input orientation. So do not use the
C        MOLECULE.INP file from this utility, neither the coordinates in the
C        end of DALTON.HES.
C        Use either the standard-orientation from Gaussian outputfile, or
C        use the converted checkpointfile (created with the command
C        "formchk checkpointfilename convertedfilename") as an input for
C        this program. Anyway, this will not lead to the geometry you 
C        specified in the input, and the Hessian will be incompatible with the
C        calculation you want to do in DALTON. (If this was already carried
C        out at the input orientation ...)
C     2) The only thing to get the geometry orientation you want is to turn
C        off symmetry. This is bloody stupid, but write to Gaussian ...
C        Then, both the Hessian and the geometry in the Test.FChk file are
C        in the input orientation, and everthing is fine.
C     #######################################################################
C
      PROGRAM FCHK2HES
C
      INTEGER HESFIL,MOLFIL,INPFIL,NATOM,NELEM,MAXAT
      PARAMETER (MAXAT=99 ,INPFIL=50,HESFIL=60,MOLFIL=70)
      DOUBLE PRECISION HESINP(3*MAXAT*(3*MAXAT+1)/2+4),
     &                 HESOUT(3*MAXAT,3*MAXAT),COORD(3*MAXAT) 
      CHARACTER STRING*80,FILENA*60
      LOGICAL READCO
      INTEGER ATNUM(0:MAXAT),NBLOCK,NLINES(MAXAT),IATOM
C
C Reading command line, to see if alternative Filename is given
C and if the coordinates must be skipped
C Commandline: FChk2HES [-n] [filename]
C -n does not read and write coordinates, default filename is Test.FChk
C
      READCO = .TRUE.
      FILENA = 'Test.FChk'
      IF (IARGC().GT.0)	THEN
          CALL GETARG (1,STRING)
          IF ((STRING(1:1).EQ.'-') .AND. (STRING(2:2).EQ.'n')) THEN
              READCO = .FALSE.
              CALL GETARG (2,STRING)
C          IF (IARGC().EQ.2) THEN
C              STRING =''
          END IF
          IF (STRING.NE.'') THEN
             FILENA = STRING
          END IF
      END IF
C
      OPEN (INPFIL,FILE=FILENA,FORM='FORMATTED',
     &STATUS='OLD',ERR=999)
C
C  Reading number of atoms
C
  5   READ (INPFIL,'(A)',END=996) STRING
      IF (INDEX(STRING,'Atomic numbers').EQ.0) GOTO  5
      NATOM = 0
      DO 8,I=59,61
         IF (STRING(I:I).NE.'') THEN 
            NATOM = NATOM*10 + (ICHAR(STRING(I:I)) - 48)
         ENDIF 
  8   CONTINUE
      IF (NATOM.GT.MAXAT) THEN
          PRINT *,'Too many atoms to read'
          GOTO 1000
      ENDIF
      NELEM = 3*NATOM*(3*NATOM+1)/2
      NCOORD = 3*NATOM
C
C    Reading atomic numbers and coordinates
C
      IF (READCO) THEN
         ATNUM(0) = 0
C We are already at atomic numbers
         DO 10,I=1,NATOM,6
             READ (INPFIL,'(6I12)',END=994) (ATNUM(I+K),K=0,5)
 10      CONTINUE
 15      READ (INPFIL,'(A)',END=995) STRING
         IF (INDEX(STRING,'Current cartesian coordinates').EQ.0) GOTO 15
         DO 20,I=1,NCOORD,5
             READ (INPFIL,'(5E16.8E2)',END=995) (COORD(I+K),K=0,4)
 20      CONTINUE
      END IF
C
C   Reading Hessian
C
 30   READ (INPFIL,'(A)',END=997) STRING
      IF (INDEX(STRING,'Cartesian Force Constants').EQ.0) GOTO 30
      DO 35,I=1,NELEM,5
         READ (INPFIL,'(5E16.8E2)',END=998) (HESINP(I+K),K=0,4)
 35   CONTINUE
      CLOSE (INPFIL)
C
C  Converting diagonal matrix to colones
C
      IOUT = 1
      DO 45,I=1,3*NATOM
         DO 40,J=1,I
             HESOUT(I,J) = HESINP(IOUT)
             HESOUT(J,I) = HESINP(IOUT)
             IOUT = IOUT + 1
 40      CONTINUE
 45   CONTINUE
C
C   Writing Hessian
C
      OPEN (HESFIL,FILE='DALTON.HES',FORM='FORMATTED',
     &STATUS='NEW',ERR=899)
      WRITE (HESFIL,'(I3)') NCOORD
      WRITE (HESFIL,'(A)') ''
      DO 60,I=1,NCOORD
         DO 50,J=1,NCOORD
            WRITE (HESFIL,'(E16.8E2)') HESOUT(I,J)
 50      CONTINUE
         WRITE (HESFIL,'(A)') ''
 60   CONTINUE
C
C  Writing coordinates in DALTON.HES
C
      IF (READCO) THEN
          DO 70,I=1,NCOORD
                WRITE (HESFIL,'(3F17.10)') COORD(I)
  70      CONTINUE 
          WRITE (HESFIL,'(A)') ''
      END IF
      CLOSE (HESFIL)
C
C  Writing coordinates in MOLECULE.INP
C
      IF (READCO) THEN
          OPEN (MOLFIL,FILE='MOLECULE.INP',FORM='FORMATTED',
     &    STATUS='NEW',ERR=799)
          WRITE (MOLFIL,'(A)') 'BASIS'
          WRITE (MOLFIL,'(A)') 'STO-3G'
          WRITE (MOLFIL,'(A)') 'This is just a dummy MOLECULE.INP '//
     &           'file that can be used'
          WRITE (MOLFIL,'(A)') 'as a starting point when '//
     &           'the Hessian is read in.'
          NBLOCK = 0
          DO 75,I=1,NATOM
               IF (ATNUM(I).NE.ATNUM(I-1)) THEN
                  NBLOCK = NBLOCK + 1
                  NLINES(NBLOCK) = 0
               END IF
               NLINES(NBLOCK) = NLINES(NBLOCK) + 1
  75      CONTINUE
          WRITE (MOLFIL,'(I5,A)') NBLOCK,'  0 0'
          IATOM = 1
          DO 85,I=1,NBLOCK
             WRITE (MOLFIL,'(I9,A,I5)') ATNUM(IATOM),'.',NLINES(I)
             DO 80,J=1,NLINES(I)
                WRITE (MOLFIL,'(A,I3.3,3F17.10)') 'X',IATOM,
     &            (COORD((IATOM-1)*3+K),K=1,3)
                IATOM = IATOM + 1
  80         CONTINUE
  85      CONTINUE
C
C          DO 80,I=0,NATOM-1
C               IF (ATNUM(I+1).NE.ATNUM(I)) THEN
C                  WRITE (MOLFIL,'(I9,A,I5)') ATNUM(I),'.',NATOM
C               END IF
C               WRITE (MOLFIL,'(A,I3.3,3F17.10)') 'X',I+1,
C     &         (COORD(I*3+J),J=1,3)
C  80      CONTINUE
          CLOSE (MOLFIL)
          PRINT *,'File conversion completed succesfully !!'
          WRITE(*,'(I2,A,I4,A)')NATOM,' Atoms and ',NELEM,
     &            ' elements of the '//
     &            'Hessian are read and written in 2 files'
      ELSE
        PRINT *,'File conversion completed succesfully !!'
        WRITE(*,'(I4,A)') NELEM, ' elements of the Hessian are '//
     &     'read and written in DALTON.HES. No MOLECULE.INP file!!'
      END IF
      GOTO 1000
C    ############################Error labels
 799  PRINT *,'The file MOLECULE.INP already exists. Please move it.'
      GOTO 1000
C
 899  PRINT *,'The file DALTON.HES already exists. Please move it.' 
      GOTO 1000
C
 994  PRINT *,'Not enough atomic numbers'
      GOTO 1000
C
 995  PRINT *,'Not enough coordinates'
      GOTO 1000
C
 996  PRINT *,'Can not find the number of atoms'
      GOTO 1000
C
 997  PRINT *,'The input file does not contain the  Hessian'
      GOTO 1000
C
 998  PRINT *,'Not enough elements of the Hessian'
      GOTO 1000
C
 999  PRINT *,'Can not find the input file'
      GOTO 1000
C #########################End of the programm
C
 1000 END



