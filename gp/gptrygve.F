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
#include <single.h>
C
#ifdef REV_LOG
c===========================================================================
cRevision 1.5-1.9  2000/05/02  hjj
c(1) changed QUIT so that it calls TSTAMP, prints total cpu and elapsed time,
c    and stops with STOP 100
c(2) removed obsolete subroutine TIMDAT
c(3) allow unit 6 in GPOPEN;
c    do not test for unit 6 in GPCLOSE (and do not close unit 6)
c
crevision 1.4 2000/04/14 07:40:36 hjj
cGPIOIN: initialize priunit.h, new defaults:  /hjaaj apr 2000
c   LUPRI = 6 (to use stdout on unix)
c   LUW4  = 6 (same as LUPRI by default)
c   NINFO = 0, NWARN = 0 (zero counters moved to here)
c   LUERR = -9xxx (let GPOPEN assign this one, could be = stderr on unix)
c   IPRERR = 0 (no debug print on LUERR)
c----------------------------
crevision 1.3 2000/04/10 10:52:18 hjj
c changed test for "no more available file numbers" in GPOPEN and DAINTF
c the old test would use IUNTAB(100) and the dimension is only 99. /HJAaJ
c----------------------------
crevision 1.2 2000/04/10 hjj
c changed initialization to -9xxx from 0, because unit 0 is stderr on
c many systems.
c===========================================================================
#endif
C  /* Deck header */
      SUBROUTINE HEADER(HEAD,IN)
      CHARACTER HEAD*(*)
#include <priunit.h>
C
      LENGTH = LEN(HEAD)
      IF (IN .GE. 0) THEN
         INDENT = IN + 1
      ELSE
         INDENT = (72 - LENGTH)/2 + 1
      END IF
      WRITE (LUPRI, '(//,80A)') (' ',I=1,INDENT), HEAD
      WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), ('-',I=1,LENGTH)
      WRITE (LUPRI, '()')
      RETURN
      END
C  /* Deck timer */
      SUBROUTINE TIMER(TEXT,TIMSTR,TIMEND)
#include <implicit.h>
      CHARACTER TEXT*6
#include <priunit.h>
C
      IF (TEXT .EQ. 'START ') THEN
         CALL GETTIM(TIMSTR,DUMMY)
      ELSE
         CALL GETTIM(TIMEND,DUMMY)
#if defined (SYS_T3D)
         TIMEND = TIMEND + 1.0D0
#endif
         TIME   = TIMEND - TIMSTR
         TIMSTR = TIMEND
         WRITE (LUPRI,'()')
         CALL TIMTXT('>>> Time used in '//TEXT//' is',TIME,LUPRI)
         WRITE (LUPRI,'()')
         CALL FLSHFO(LUPRI)
      END IF
      RETURN
      END
C  /* Deck timpri */
      SUBROUTINE TIMPRI(TEXT,TIME,TIMALL)
#include <implicit.h>
      CHARACTER TEXT*6
      PARAMETER (HUN = 100.0D00)
#include <priunit.h>
      SAVE TREST
      DATA TREST /0.0D0/
      IF (TEXT .EQ. 'REST  ') THEN
         TIME = TREST
C 980826-hjaaj: initialize TREST for next abacus call
         TREST = 0.0D0
      END IF
      ITIME = NINT(HUN*TIME/TIMALL)
      IF (ITIME .GT. 0 .AND. INT(TIME) .GT. 0) THEN
         MINUTE = INT(TIME)/60
         IHOURS = MINUTE/60
         MINUTE = MINUTE - 60*IHOURS
         ISECND = NINT(TIME) - 3600*IHOURS - 60*MINUTE
         WRITE(LUPRI,100) TEXT, IHOURS, MINUTE, ISECND, ITIME
      ELSE
         TREST = TREST + TIME
      END IF
  100 FORMAT(1X,A6,'     ',I2.2,':',I2.2,':',I2.2,5X,I3,' %')
      RETURN
      END
C  /* Deck titler */
      SUBROUTINE TITLER(HEAD,A,IN)
      CHARACTER HEAD*(*), A
#include <priunit.h>
C
      LENGTH = LEN(HEAD)
      IF (IN .EQ. 200) THEN
         LENGTH = LENGTH + 2
      ELSE IF (IN .GE. 100) THEN
         MARG = IN - 100
         IF (MARG .GT. 0) MARG = MARG + 1
         LENGTH = LENGTH + 2*MARG
      END IF
      IF (IN .GE. 0 .AND. IN .LT. 100) THEN
         INDENT = IN + 1
      ELSE
         INDENT = MAX(1,(72 - LENGTH)/2 + 1)
      END IF
      IF (IN .EQ. 200) THEN
         WRITE (LUPRI, '(//,80A)')
     *      (' ',I=1,INDENT), '.', ('-',I=1,LENGTH), '.'
         WRITE (LUPRI,'(80A)') (' ',I=1,INDENT),'| ', HEAD, ' |'
         WRITE (LUPRI, '(80A)')
     *      (' ',I=1,INDENT), '`', ('-',I=1,LENGTH), ''''
      ELSE IF (IN .EQ. 100) THEN
         WRITE (LUPRI, '(//,80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), HEAD
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      ELSE IF (IN .GT. 100) THEN
         WRITE (LUPRI, '(//,80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT),
     *      (A,I=1,MARG-1), ' ', HEAD, ' ', (A,I=1,MARG-1)
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      ELSE
         WRITE (LUPRI, '(//,80A)') (' ',I=1,INDENT), HEAD
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      END IF
      WRITE (LUPRI, '()')
      CALL FLSHFO(LUPRI)
      RETURN
      END
C  /* Deck around */
      SUBROUTINE AROUND(HEAD)
      CHARACTER HEAD*(*)
#include <priunit.h>
      LNG = LEN(HEAD) + 2
      IND = (72 - LNG)/2 + 1
      WRITE (LUPRI, '(//,80A)') (' ',I=1,IND), '+', ('-',I=1,LNG), '+'
      WRITE (LUPRI, '(80A)')    (' ',I=1,IND), '! ', HEAD, ' !'
      WRITE (LUPRI, '(80A)')    (' ',I=1,IND), '+', ('-',I=1,LNG), '+'
C     WRITE (LUPRI, '(//,80A)') (' ',I=1,IND), '.', ('-',I=1,LNG), '.'
C     WRITE (LUPRI, '(80A)')    (' ',I=1,IND), '| ', HEAD, ' |'
C     WRITE (LUPRI, '(80A)')    (' ',I=1,IND), '`', ('-',I=1,LNG), ''''
      WRITE (LUPRI, '()')
      RETURN
      END
C  /* Deck stopit */
      SUBROUTINE STOPIT(SUB,PLACE,INT1,INT2)
#include <implicit.h>
#include <priunit.h>
      CHARACTER*(*) SUB, PLACE
      WRITE (LUPRI,'(//3A)') ' Work space exceeded in subroutine ',
     *                         SUB,'.'
      IF ((LEN(PLACE) .GT. 1) .OR. (PLACE .NE. ' ')) THEN
         WRITE (LUPRI,'(/2A)') ' Location: ',PLACE
      END IF
      JNT1  = ABS(INT1)
      JNT2  = ABS(INT2)
      LWRKR = MAX(JNT1,JNT2)
      LWRKA = MIN(JNT1,JNT2)
      IF (INT1 .LT. 0 .OR. INT2 .LT. 0) THEN
         WRITE (LUPRI,'(/A,I10)  ') ' Space required  >',LWRKR
      ELSE
         WRITE (LUPRI,'(/A,I10)  ') ' Space required  =',LWRKR
      END IF
      WRITE (LUPRI,'( A,I10)') ' Space available =',LWRKA
      CALL QUIT('Work space exceeded.')
      END
C  /* Deck quit */
      SUBROUTINE QUIT(TEXT)
#include <implicit.h>
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
#include <priunit.h>
      CHARACTER TEXT*(*)
C
C     Stamp date and time and hostname to output
C     (if output unit has been defined yet)
C
      IF (LUPRI .GE. 0) THEN
         LUPRIN = LUPRI
      ELSE IF (LUERR .GE. 0) THEN
         LUPRIN = LUERR
      ELSE
         LUPRIN = 6
      END IF
      CALL TSTAMP('  --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---',
     &   LUPRIN)
      WRITE (LUPRIN,'(/2A/)') ' Reason: ',TEXT
#if defined (SYS_CRAY) || defined (SYS_UNIX) || defined (SYS_AIX) || defined (SYS_HPUX) || defined (SYS_IRIX) || defined (SYS_NEC) || defined (SYS_DEC) || defined (SYS_T3D) || defined (SYS_LINUX) || defined (SYS_SUN) || defined (SYS_HPUX) || defined (SYS_HAL) || defined (SYS_T90)
C     Write to stderr
      WRITE (0,'(/A/1X,A)')
     *   '  --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---',TEXT
#endif
      CALL GETTIM(CTOT,WTOT)
      CALL TIMTXT('>>>> Total CPU  time used in DALTON:',CTOT,LUPRIN)
      CALL TIMTXT('>>>> Total wall time used in DALTON:',WTOT,LUPRIN)
      CALL QTRACE(LUPRIN)
      CALL FLSHFO(LUPRIN)
C
#if defined (VAR_MPI)
      IERR = 100
C     ... define exit code to the environment for MPI_ABORT
      CALL MPI_ABORT(MPI_COMM_WORLD,IERR,IERROR)
#endif
#if defined (VAR_PVM)
      CALL PVMXHALT
#endif
C
      CALL TRACE
      STOP 100
      END
C  /* Deck opendx */
      SUBROUTINE OPENDX (LUDX,NAME,NELEM,STATUS,OLDDX)
C
C 15-Jun-1985 hjaaj
C
C Revisions :  9-Dec-1987 hjaaj (Alliant version)
C
C Purpose:
C   Open files for direct access through WRITDX and READDX routines.
C   The ....DX routines enables direct access, even when the number
C   of elements per record (the logical record length) is greater
C   than the maximum physical record length.
C   >>> THIS IS MACHINE DEPENDENT <<<
C
C Input:
C  LUDX     file unit number
C  NELEM    number of integer words per logical record
C  STATUS   file status: 'OLD', 'NEW', or 'UNKNOWN'
C
C Output:
C  OLDDX    logical, true if old LUDX file was opened
C
C
      CHARACTER*(*) NAME, STATUS
      LOGICAL OLDDX, FEXIST
#include <priunit.h>
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
C
C     CRAY has 8 byte integers.
C
      LRECL  = 8*NELEM
C
      IF (STATUS .EQ. 'NEW') THEN
         INQUIRE(FILE=NAME,EXIST=FEXIST)
         IF (FEXIST) THEN
            OPEN(UNIT=LUDX,STATUS='OLD',FILE=NAME)
            CLOSE(UNIT=LUDX,STATUS='DELETE')
         END IF
         GO TO 300
      END IF
C
      IF (STATUS .NE. 'OLD' .AND. STATUS .NE. 'UNKNOWN') GO TO 9000
C
C     OPEN OLD FILE
C
         OPEN(LUDX,FILE=NAME,STATUS='OLD',FORM='UNFORMATTED',ERR=300,
     *        ACCESS='DIRECT',RECL=LRECL)
         OLDDX = .TRUE.
      GO TO 600
C
  300 CONTINUE
      IF (STATUS .EQ. 'OLD') GO TO 9100
C
C     OPEN NEW FILE
C
         OPEN(LUDX,FILE=NAME,STATUS='NEW',FORM='UNFORMATTED',
     *        ACCESS='DIRECT',ERR=301,IOSTAT=IOS,RECL=LRECL)
         OLDDX = .FALSE.
C
C     ERROR HANDLING (JAA 1997-06-06)
C
 301     IF (IOS .GT. 0) THEN
            WRITE(102,*) '--> ERROR (OPENDX) NEW FILE OPEN FAILED'
            WRITE(102,*) '    FILE=',NAME,', UNIT=',LUDX,
     $           ', IOSTAT=',IOS
            IF (IOS .LT. 1000) THEN
               WRITE(102,'(A,I3.3,A)') '    (try explain sys-',IOS,')'
            ELSE IF (IOS .LT. 6000) THEN
               WRITE(102,'(A,I4,A)') '    (try explain lib-',IOS,')'
            ENDIF
            CALL FLUSH(102,IOS)
            CALL TRBK()
#if defined(SYS_T3D)
            CALL STOP_ALL()
#else
            CALL QUIT('OPENDX: new file open failed')
#endif
         ENDIF
  600 CONTINUE
#endif
#if defined (SYS_AIX) || defined (SYS_HPUX) || defined (SYS_DEC) || defined (SYS_IRIX) || defined (SYS_NEC) || defined (SYS_LINUX) || defined (SYS_SUN) || defined (SYS_HAL)
C
C     Assume same for IBM-AIX--  1-Oct-1990 hjaaj
C     Assume same for HP-UX  -- 21-Aug-1991 hjaaj
C     Assume same for DEC    -- 21-May-1992 hjaaj
C     Assume same for IRIX   --  3-Feb-1994 hjaaj; except recl in *4 words
C     Assume same for PARAGON-- 13-Oct-1994 hjaaj
C
#if (defined (SYS_IRIX) && defined(NO_F90)) || defined (SYS_NEC)
      LRECL  = NELEM
#else
#ifdef VAR_INT64
      LRECL  = 8*NELEM
#else 
      LRECL  = 4*NELEM
#endif
#endif
C
      IF (STATUS .EQ. 'NEW') THEN
C        some systems abort if we try to open as NEW a file,
C        which already exist
         INQUIRE(FILE=NAME,EXIST=FEXIST)
         IF (FEXIST) THEN
            OPEN(UNIT=LUDX,STATUS='OLD',FILE=NAME)
            CLOSE(UNIT=LUDX,STATUS='DELETE')
         END IF
         GO TO 300
      END IF
C
      IF (STATUS .NE. 'OLD' .AND. STATUS .NE. 'UNKNOWN') GO TO 9000
C
C     OPEN OLD FILE
C
         OPEN(LUDX,FILE=NAME,STATUS='OLD',FORM='UNFORMATTED',ERR=300,
     *        ACCESS='DIRECT',RECL=LRECL)
         OLDDX = .TRUE.
      GO TO 600
C
  300 CONTINUE
      IF (STATUS .EQ. 'OLD') GO TO 9100
C
C     OPEN NEW FILE
C
         OPEN(LUDX,FILE=NAME,STATUS='NEW',FORM='UNFORMATTED',
     *        ACCESS='DIRECT',RECL=LRECL)
         OLDDX = .FALSE.
  600  CONTINUE
#endif
#if !defined (SYS_CRAY) && !defined (SYS_AIX) && !defined (SYS_HPUX) && !defined (SYS_DEC) && !defined (SYS_IRIX) && !defined (SYS_NEC) && !defined (SYS_T3D) && !defined (SYS_LINUX) && !defined (SYS_SUN) && !defined (SYS_HAL) && !defined (SYS_T90)
 >>>>> insert appropriate OPEN statements in OPENDX.
#endif
      RETURN
C
C error branches
C
 9000 CONTINUE
      WRITE (LUPRI,'(//2A/2A/A,I5)')
     *   ' --> ERROR (OPENDX) INVALID STATUS KEYWORD: ',STATUS,
     *   '                    FILE NAME   =',NAME,
     *   '                    FILE NUMBER =',LUDX
      CALL QTRACE(LUPRI)
      CALL QUIT('ERROR (OPENDX) INVALID STATUS KEYWORD')
C
 9100 CONTINUE
      WRITE (LUPRI,'(//A/2A/A,I5/A)')
     *   ' --> ERROR (OPENDX) OLD FILE NOT FOUND',
     *   '                    FILE NAME   =',NAME,
     *   '                    FILE NUMBER =',LUDX,
     *   ' --- or wrong record length on old file.'
      CALL QTRACE(LUPRI)
      CALL QUIT('ERROR (OPENDX) FILE NOT FOUND')
C
C end of OPENDX
C
      END
C  /* Deck finddx */
      LOGICAL FUNCTION FINDDX(LU,I,LEN,IVEC)
C
C 27-Jun-1985 Hans Jorgen Aa. Jensen
C
C For direct access find record,
C when LEN may be greater than maximum record length.
C
      INTEGER IVEC(LEN)
      READ (LU, REC=I, IOSTAT=IOS, ERR=900) IVEC
      IF (IOS .NE. 0) GO TO 900
      FINDDX = .TRUE.
      RETURN
C
  900 CONTINUE
      FINDDX = .FALSE.
      RETURN
      END
C  /* Deck readdx */
      SUBROUTINE READDX(LU,I,LEN,IVEC)
C
C 30-Apr-1985 Hans Jorgen Aa. Jensen
C
C For direct access read
C
#include <implicit.h>
#include <priunit.h>
#include <2gbdef.h>
#include <dummy.h>
      LOGICAL OLDDX
      CHARACTER*80 FNNAME, FNNM2
      INTEGER IVEC(LEN)
#include <inftap.h>
#include <chrnos.h>
C
      IF (LEN .EQ. 0) RETURN
#if defined (VAR_SPLITFILES)
      INQUIRE(UNIT=LU,RECL=LRECL)
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
      LRECL = LRECL/8
#else
#if !defined (SYS_NEC) && !defined (SYS_IRIX)
#ifdef VAR_INT64
      LRECL = LRECL/8
#else
      LRECL = LRECL/4
#endif
#endif
#endif
      IRECPF = I2GB/LRECL
C
C     Error branch for very large record lengths
C
      IF (IRECPF .EQ. 0) THEN
         WRITE (LUPRI,'(/A,I14,A)') 
     &        ' ERROR: The requested record length of ', LRECL, 
     &        ' is larger than 2Gb'
         WRITE (LUPRI,'(/A)') ' The program cannot write even a '//
     &        'single record to file!'
         CALL QUIT('Record length of file too long for this computer')
      END IF
C
      IF (I .GT. IRECPF) THEN
C
C     Ooops, this record is not in this file
C
         LFILNM = I/IRECPF - 1
         IF (LFILNM .GT. 9) THEN
            WRITE (LUPRI,'(/A)') ' DALTON needs to read from a file '//
     &           ' split more than 11 times.',
     &           ' This is currently not supported'
            CALL QUIT('Too many splittings of a file')
         END IF
         IRECT = MOD(I,IRECPF)
         IF (IRECT .EQ. 0) IRECT = IRECPF
         INQUIRE(UNIT=LU,NAME=FNNAME)
         LN = 1
 10      CONTINUE
         IF (FNNAME(LN:LN) .NE. ' ') THEN
            LN = LN + 1
            GOTO 10
         END IF
         LN = LN - 1
         LUBKP = LU
         CALL GPCLOSE(LU,'KEEP')
         LU = LUBKP
         FNNM2 = FNNAME(1:LN)//'-'//CHRNOS(LFILNM)
         CALL GPOPEN(LU,FNNM2(1:(LN+2)),'UNKNOWN','DIRECT',' ',LRECL,
     &               OLDDX)
      ELSE
         IRECT = I
      END IF
#else
      IRECT = I
#endif
      READ (LU, REC = IRECT) IVEC
#if defined (VAR_SPLITFILES)
      IF (I .GT. IRECPF) THEN
         LUBKP = LU
         CALL GPCLOSE(LU,'KEEP')
         LU = LUBKP
         CALL GPOPEN(LU,FNNAME(1:LN),'UNKNOWN','DIRECT',' ',LRECL,OLDDX)
      END IF
#endif
      RETURN
      END
C  /* Deck writdx */
      SUBROUTINE WRITDX(LU,I,LEN,IVEC)
C
C 30-Apr-1985 Hans Jorgen Aa. Jensen
C
#include <implicit.h>
#include <dummy.h>
#include <2gbdef.h>
#include <priunit.h>
      LOGICAL OLDDX
      CHARACTER*80 FNNAME, FNNM2
      INTEGER IVEC(LEN)
#include <inftap.h>
#include <chrnos.h>
C
#if defined (VAR_SPLITFILES)
      INQUIRE(UNIT=LU,RECL=LRECL)
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
      LRECL = LRECL/8
#else
#if !defined (SYS_NEC) && !defined (SYS_IRIX)
#ifdef VAR_INT64
      LRECL = LRECL/8
#else
      LRECL = LRECL/4
#endif
#endif
#endif
      IRECPF = I2GB/LRECL
C
C     Error branch for very large record lengths
C
      IF (IRECPF .EQ. 0) THEN
         WRITE (LUPRI,'(/A,I14,A)') 
     &        ' ERROR: The requested record length of ', LRECL, 
     &        ' is larger than 2Gb'
         WRITE (LUPRI,'(/A)') ' The program cannot write even a '//
     &        'single record to file!'
         CALL QUIT('Record length of file too long for this computer')
      END IF
C
      IF (I .GT. IRECPF) THEN
C
C     Ooops, this record will not fit in this file
C
         LFILNM = I/IRECPF - 1
         IF (LFILNM .GT. 9) THEN
            WRITE (LUPRI,'(/A)') ' DALTON need to split a file more '//
     &           ' than 11 times.', ' This is currently not supported'
            CALL QUIT('Too many splittings of a file')
         END IF
         IRECT = MOD(I,IRECPF)
         IF (IRECT .EQ. 0) IRECT = IRECPF
         INQUIRE(UNIT=LU,NAME=FNNAME)
         LN = 1
 10      CONTINUE
         IF (FNNAME(LN:LN) .NE. ' ') THEN
            LN = LN + 1
            GOTO 10
         END IF
         LN = LN - 1
         LUBKP = LU
         CALL GPCLOSE(LU,'KEEP')
         LU = LUBKP
         FNNM2 = FNNAME(1:LN)//'-'//CHRNOS(LFILNM)
         LN2 = LN + 2
         CALL GPOPEN(LU,FNNM2(1:LN2),'UNKNOWN','DIRECT',' ',LRECL,
     &               OLDDX)
      ELSE
         IRECT = I   
      END IF
#else
      IRECT = I
#endif
      WRITE (LU, REC = IRECT) IVEC
#if defined (VAR_SPLITFILES)
      IF (I .GT. IRECPF) THEN
         LUBKP = LU
         CALL GPCLOSE(LU,'KEEP')
         LU = LUBKP
         CALL GPOPEN(LU,FNNAME(1:LN),'UNKNOWN','DIRECT',' ',LRECL,OLDDX)
      END IF
#endif
      RETURN
      END
C  /* Deck aolab4 */
      SUBROUTINE AOLAB4(IINDPK,NMAX,NIBUF,NBITS,IINDX4,N)
C
C     Written by Henrik Koch 22-Nov-1991 (as AOLABE)
C     Generalized 25-Oct-1993 hjaaj
C                 31-Jul-1996 hjaaj (new NIBUF parameter)
C
C     Unpack 4 integer indices packed in NIBUF *4 integer(s).
C
#include <implicit.h>
      DIMENSION IINDPK(*), IINDX4(4,*)
#include <ibtdef.h>
C
      N = IINDPK(NIBUF*NMAX+1)
C
      IF (NIBUF .EQ. 1) THEN
      IF (NBITS .EQ. 8) THEN
         DO 100 I = 1,N
            LABEL = IINDPK(I)
            IINDX4(1,I) = IBTAND(IBTSHR(LABEL,24),IBT08)
            IINDX4(2,I) = IBTAND(IBTSHR(LABEL,16),IBT08)
            IINDX4(3,I) = IBTAND(IBTSHR(LABEL, 8),IBT08)
            IINDX4(4,I) = IBTAND(       LABEL,    IBT08)
  100    CONTINUE
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
      ELSE IF (NBITS .EQ. 16) THEN
         DO 200 I = 1,N
            LABEL = IINDPK(I)
            IINDX4(1,I) = IBTAND(IBTSHR(LABEL,48),IBT16)
            IINDX4(2,I) = IBTAND(IBTSHR(LABEL,32),IBT16)
            IINDX4(3,I) = IBTAND(IBTSHR(LABEL,16),IBT16)
            IINDX4(4,I) = IBTAND(       LABEL,    IBT16)
  200    CONTINUE
#endif
      ELSE
         CALL QUIT('AOLAB4 error')
      END IF
      ELSE
         DO 300 I = 1,N
            LABEL = IINDPK(I)
            IINDX4(1,I) = IBTAND(IBTSHR(LABEL,16),IBT16)
            IINDX4(2,I) = IBTAND(       LABEL    ,IBT16)
            LABEL = IINDPK(NMAX+I)
            IINDX4(3,I) = IBTAND(IBTSHR(LABEL,16),IBT16)
            IINDX4(4,I) = IBTAND(       LABEL    ,IBT16)
  300    CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck GPOPEN */
      SUBROUTINE GPOPEN(LUNIT,FILEIN,STATIN,ACCEIN,FORMIN,NELEM,OLDDX)
C
C 23-Feb-2000 K.Ruud
C
C Purpose:
C   General purpose routine for opening files in the Dalton program.
C   The routine will dynamically allocate unit numbers that will
C   become available again when the file is closed using GPCLOSE
C
C   Direct access files are provided with a UNIT-number, but the rest
C   of the file opening process is taken care of by a call to OPENDX.
C
C   These files are strongly machine dependent, although care has been 
C   taken to avoid using unit numbers that are illegal or reserved on a
C   particular architecture. However, the use of GPOPEN and GPCLSE should
C   remove much of the machine dependence in the rest of the Dalton 
C   program
C
C Input:
C  LUNIT    Suggested unit number (OPTIONAL, but is mandatory if an unnamed
C           file is reopened after having been closed with STATUS='KEEP')
C           Otherwise it is recommended to not assign this.
C  FILEIN   Suggested name for the file (OPTIONAL, but strongly recommended)
C  STATIN   Suggested status of the file (OPTIONAL and maybe not recommended)
C  ACCEIN   Access method. 'DIRECT' or 'SEQUENTIAL'. Default is sequential.
C  FORMIN   Formatted or unformatted file format. Default is 'UNFORMATTED'
C
C Direct access-specific input:
C
C  NELEM    number of integer words per logical record
C  NREC     number of logical records
C
C
C Output:
C  LUNIT    Assigned file unit number
C
C  Direct access-specific output (through 
C  LRDX     physical record length (in integers)
C  NBDX     number of physical records per logical record
C  OLDDX    logical, true if old LUDX file was opened
C
C
#include <implicit.h>
#include <priunit.h>
#include <chrnos.h>
      CHARACTER*80 FILENM, STATUF, ACCESF, FORMF
      CHARACTER*(*) FILEIN, STATIN, ACCEIN, FORMIN
      CHARACTER OUTFIL*20, WRKDIR*60
      INTEGER FILELN, STATLN, ACCELN, FORMLN
      LOGICAL OLDDX, FEXST, FIRST
      COMMON /UNITAR/ IUNTAB(1:99)
#include <gnrinf.h>
c maxorb for infpar.h
#include <maxorb.h>
#include <infpar.h>
      DATA FIRST /.TRUE./
      SAVE FIRST
C
C     Initialization
C     
      IF (FIRST) THEN
         CALL IZERO(IUNTAB,99)
#if defined (SYS_HPUX)
C
C     This avoids us writing to unit 7 of HP-systems (which corresponds to
C     standard error on this machine)
C
         IUNTAB(7) = 1
#endif
         FIRST = .FALSE.
      END IF
C
C     We first deal with the unit number
C
      FILELN = LEN(FILEIN)
      STATLN = LEN(STATIN)
      ACCELN = LEN(ACCEIN)
      FORMLN = LEN(FORMIN)
      FILENM(1:FILELN) = FILEIN(1:FILELN)
      STATUF(1:STATLN) = STATIN(1:STATLN)
      ACCESF(1:ACCELN) = ACCEIN(1:ACCELN)
      FORMF (1:FORMLN) = FORMIN(1:FORMLN)
      IF ((LUNIT .LT. 1) .OR. (LUNIT .GT. 99)) THEN
C
C     Unit number left unassigned, we get to decide!
C
         IUN = 0
 10      CONTINUE
         IUN = IUN + 1
         IF (IUN .GT. 99) GOTO 9001
         IF ((IUN .EQ. 5) .OR. (IUN .EQ. 6))
     &        GOTO 10
         IF (IUNTAB(IUN) .NE. 0) GOTO 10
         LUNIT = IUN
      ELSE
C
C     The user has requested a specific unit number. We don't 
C     quite trust the user, so we check that
C     1) It is not unit 5, or 6
C     2) The file either has been closed with status='KEEP' or is
C        not currently in use
C
         IF ((LUNIT .EQ. 5) .OR. (LUNIT .EQ. 6)) GOTO 9002
         IF (IUNTAB(LUNIT) .EQ. 1) GOTO 9003
      END IF
#if defined (VAR_MPI)
      IF (MYNUM .EQ. 0) THEN
#endif
      IF (FILENM(1:FILELN) .EQ. 'DALTON.OUT') THEN
         OUTFIL = '                    '
#if defined (SYS_T3D) || defined (SYS_T90)
         CALL PXFGETENV ('OUTFIL',6,OUTFIL,LOUTFL,IERR)
#else 
         CALL GETENV ('OUTFIL',OUTFIL)
#endif
         IF (OUTFIL(1:1) .NE. ' ') THEN
            LENOUT = 0
            DO I = 1, 20
               IF (OUTFIL(I:I) .EQ. ' ') GO TO 43
               LENOUT = LENOUT + 1
            END DO
 43         CONTINUE
#if defined (SYS_T3D) || defined (SYS_T90)
            CALL PXFGETENV ('WRKDIR',6,WRKDIR,LWRKDR,IERR)
#else 
            CALL GETENV ('WRKDIR',WRKDIR)
#endif
            LENWRK = 0
            DO I = 1, 60
               IF (WRKDIR(I:I) .EQ. ' ') GO TO 44
               LENWRK = LENWRK + 1
            END DO
 44         CONTINUE
            FILELN = LENWRK + LENOUT + 1
            FILENM(1:FILELN) = WRKDIR(1:LENWRK)//'/'//
     &                         OUTFIL(1:LENOUT)
         END IF
      END IF
#if defined (VAR_MPI)
      END IF
#endif
      IF (STATLN .EQ. 1) THEN
         STATLN = 7
         STATUF(1:7) = 'UNKNOWN'
      END IF
C
C     We've got a file number now, get rid of the direct access case
C
      IF (ACCESF(1:6) .EQ. 'DIRECT') THEN
         IF (FILELN .EQ. 1) THEN
#if defined (VAR_MPI)
            IF (MYNUM .EQ. 0) THEN
#endif
               FILELN = 9
               FILENM(1:9) = 'UNIT'//CHRNOS(LUNIT/10)//
     &                        CHRNOS(MOD(LUNIT,10))//'.DA'
#if defined (VAR_MPI)
            ELSE
               FILELN = 13
               FILENM(1:9) = 'UNIT'//CHRNOS(LUNIT/10)//
     &                        CHRNOS(MOD(LUNIT,10))//'.DA'//
     &                        '#'//CHRNOS(MYNUM/100)//
     &                        CHRNOS(MOD(MYNUM,100)/10)//
     &                        CHRNOS(MOD(MOD(MYNUM,100),10))
            END IF
         ELSE
            IF (FILENM(FILELN-3:FILELN-3) .NE. '#' .AND.
     &          MYNUM .GT. 0) THEN
               FILENM(FILELN+1:FILELN+4) = '#'//CHRNOS(MYNUM/100)//
     &                        CHRNOS(MOD(MYNUM,100)/10)//
     &                        CHRNOS(MOD(MOD(MYNUM,100),10))
               FILELN = FILELN + 4
            END IF
#endif
         END IF
         IF (IPRERR.GT.1) WRITE (LUERR,*) 'GPOPEN calling OPENDX ',
     &      LUNIT,' ',FILENM(1:FILELN),NELEM,' ',STATUF(1:STATLN)
         CALL OPENDX(LUNIT,FILENM(1:FILELN),NELEM,STATUF(1:STATLN),
     &               OLDDX)
      ELSE
         IF (STATUF(1:3) .EQ. 'OLD' .AND. IUNTAB(LUNIT) .EQ. 0) THEN
C
C     This better be a file with a name, and it better exist
C
            IF (FILELN .EQ. 1) GOTO 9005
#if defined (VAR_MPI)
            IF (FILENM(FILELN-3:FILELN-3) .NE. '#' .AND.
     &           MYNUM .GT. 0 .AND. FILENM(1:3) .NE. BASDIR(1:3) 
     &           .AND. FILENM(FILELN-5:FILELN) .NE. 'NUMCAL') THEN
C
C     The last test will not be necessary anymore when we get
C     locking RMA operations in MPI2
C
               FILENM(FILELN+1:FILELN+4) = '#'//CHRNOS(MYNUM/100)//
     &                        CHRNOS(MOD(MYNUM,100)/10)//
     &                        CHRNOS(MOD(MOD(MYNUM,100),10))
               FILELN = FILELN + 4
            END IF
#endif
            INQUIRE(FILE=FILENM(1:FILELN),EXIST=FEXST,IOSTAT=IOS)
            IF (.NOT. FEXST) GOTO 9007
            IF (FORMLN .EQ. 1) THEN
               FORMLN = 11
               FORMF(1:FORMLN)='UNFORMATTED'
            END IF
            IF (IPRERR.gt.1) write (LUERR,*) 'GPOPEN seq OLD ',LUNIT,
     &         ' ',FILENM(1:FILELN),' ',FORMF(1:FORMLN)
            OPEN(UNIT=LUNIT,FILE=FILENM(1:FILELN),STATUS='OLD',
     &           ACCESS='SEQUENTIAL',FORM=FORMF(1:FORMLN),ERR=9004,
     &           IOSTAT=IOS)
         ELSE
            IF (STATUF(1:3) .EQ. 'NEW' .AND. IUNTAB(LUNIT) .NE. 0)
     &           GOTO 9006
            IF (FILELN .EQ. 1) THEN
#if defined (VAR_MPI)
               IF (MYNUM .EQ. 0) THEN
#endif
                  FILELN = 6
                  FILENM = 'UNIT'//CHRNOS(LUNIT/10)//
     &                     CHRNOS(MOD(LUNIT,10))
#if defined (VAR_MPI)
               ELSE
                  FILELN = 10
                  FILENM = 'UNIT'//CHRNOS(LUNIT/10)//
     &                      CHRNOS(MOD(LUNIT,10))//
     &                      '#'//CHRNOS(MYNUM/100)//
     &                      CHRNOS(MOD(MYNUM,100)/10)//
     &                      CHRNOS(MOD(MOD(MYNUM,100),10))
               END IF
            ELSE
               IF (FILENM(FILELN-3:FILELN-3) .NE. '#' .AND.
     &             MYNUM .GT. 0) THEN
                  FILENM(FILELN+1:FILELN+4) = '#'//CHRNOS(MYNUM/100)//
     &                                     CHRNOS(MOD(MYNUM,100)/10)//
     &                                  CHRNOS(MOD(MOD(MYNUM,100),10))
                  FILELN = FILELN + 4
               END IF
#endif
            END IF
            IF (FORMLN .EQ. 1) THEN
               FORMLN = 11
               FORMF(1:11)='UNFORMATTED'
            END IF
            IF (STATUF(1:7) .EQ. 'SCRATCH') THEN
               IF (IPRERR.gt.1) write (LUERR,*) 'GPOPEN seq scratch ',
     &            LUNIT,' ',FORMF(1:FORMLN)
               OPEN(UNIT=LUNIT,STATUS='SCRATCH',
     &              ACCESS='SEQUENTIAL',FORM=FORMF(1:FORMLN),
     &              ERR=9004,IOSTAT=IOS)
            ELSE
               IF (IPRERR.gt.1) write (LUERR,*) 'GPOPEN seq ',
     &            STATUF(1:STATLN),LUNIT,' ',FILENM(1:FILELN),
     &            ' ',FORMF(1:FORMLN)
               IF(STATUF(1:STATLN) .EQ. 'NEW') THEN
Chjaaj/may2000:... if filename already exists the program will abort
C                  thus we inquire first and delete if necessary.
C                  This will often be the case if we restart a
C                  calculation.
                  INQUIRE(FILE=FILENM(1:FILELN),EXIST=FEXST,IOSTAT=IOS)
                  IF (FEXST) THEN
                     IF (IPRERR.gt.1) write (LUERR,*) 'GPOPEN INFO: '//
     &               'this file exists and is deleted first.'
                     OPEN(UNIT=LUNIT,FILE=FILENM(1:FILELN),
     &                   STATUS='OLD',ACCESS='SEQUENTIAL',
     &                   FORM=FORMF(1:FORMLN),ERR=9004,IOSTAT=IOS)
                     CLOSE(UNIT=LUNIT,STATUS='DELETE')
                  END IF
               END IF
               OPEN(UNIT=LUNIT,FILE=FILENM(1:FILELN),
     &              STATUS=STATUF(1:STATLN),ACCESS='SEQUENTIAL',
     &              FORM=FORMF(1:FORMLN),ERR=9004,IOSTAT=IOS)
            END IF
         END IF
      END IF
  100 IUNTAB(LUNIT) = 1
      RETURN
C     
C error branches
C
 9001 CONTINUE
      WRITE (6,'(//A/A/A//A)')
     &   '--> ERROR (GPOPEN) NO MORE AVAILABLE FILENUMBERS!',
     &   '--> THIS CALCULATION EITHER NEEDS TOO MANY SIMULTANEOUS '//
     &   'FILES OR', 
     &   '--> SOMEBODY HAS FORGOTTEN TO CLOSE FILES IN THE SOURCE CODE',
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
      CALL QTRACE(6)
      CALL QUIT('ERROR (GPOPEN) NO MORE FILE NUMBERS')
C
 9002 CONTINUE
      WRITE (6,'(//A/A,I3/A//2A//A)')
     &   '--> ERROR (GPOPEN) TRYING TO OPEN AN ILLEGAL FILE NUMBER',
     &   '--> SOMEBODY HAS TRIED TO OPEN UNIT',LUNIT,
     &   '--> THE PROGRAM DOES NOT ALLOW THE USE OF THIS RESERVED'//
     &   'UNIT NUMBER',
     &   '--> Name of offending file (if any): ',FILENM(1:FILELN),
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
      CALL QTRACE(6)
      CALL QUIT('ERROR (GPOPEN) ILLEGAL FILE NUMBER REQUESTED')
C
 9003 CONTINUE
      WRITE (6,'(//A/A//A,I5/2A//A)')
     &   '--> ERROR (GPOPEN) TRYING TO OPEN A FILE ALREADY IN USE',
     &   '--> SOMEBODY IS TRYING TO USE A FILENUMBER THAT IS '//
     &   'ALREADY IN USE',
     &   '--> Offending UNIT number: ',LUNIT,
     &   '--> Name of offending file (if any): ', FILENM(1:FILELN),
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
      CALL QTRACE(6)
      CALL QUIT('ERROR (GPOPEN) TRYING TO OPEN A FILE ALREADY IN USE')
C
 9004 CONTINUE
      WRITE (6,'(//A,I3/2A/A,I7)')
     &   '--> ERROR (GPOPEN) UPON TRYING TO OPEN FILE ON UNIT',LUNIT,
     &   '--> with filename ',FILENM(1:FILELN),
     &   '--> IOSTAT ERROR CODE RETURNED ',IOS
      CALL QTRACE(6)
      CALL QUIT('ERROR (GPOPEN) UPON OPENING A FILE')
C
 9005 CONTINUE
      WRITE (6,'(//A/A/A//A)')
     &   '--> ERROR (GPOPEN) TRYING TO OPEN A NON-EXISTING OLD FILE',
     &   '--> A FILE HAS BEEN SPECIFIED TO BE OLD, BUT THE TABLE',
     &   '--> ENTRY INDICATES THAT IT DOES NOT EXIST',
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
      CALL QTRACE(6)
      CALL QUIT('ERROR (GPOPEN) TRYING TO OPEN A NON-EXISTING '//
     &          'FILE AS OLD')
C
 9006 CONTINUE
      WRITE (6,'(//A/A/A//A,I5/2A//A)')
     &   '--> ERROR (GPOPEN) TRYING TO OPEN AN EXISTING NEW FILE',
     &   '--> A FILE HAS BEEN SPECIFIED TO BE NEW, BUT THE TABLE',
     &   '--> ENTRY INDICATES THAT IT ALREADY EXISTS',
     &   '--> Offending UNIT number: ',LUNIT,
     &   '--> Name of offending file (if any): ',FILENM(1:FILELN),
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
      CALL QTRACE(6)
      CALL QUIT('ERROR (GPOPEN) TRYING TO OPEN AN EXISTING '//
     &          'FILE AS NEW')
C
 9007 CONTINUE
      WRITE (6,'(//A/A/A//A,I5/2A//A)')
     &   '--> ERROR (GPOPEN) TRYING TO OPEN A NON-EXISTING OLD FILE',
     &   '--> A FILE HAS BEEN SPECIFIED TO BE OLD,',
     &   '--> BUT THE FILE DOES NOT EXIST',
     &   '--> Offending UNIT number: ',LUNIT,
     &   '--> Name of offending file (if any): ',FILENM(1:FILELN),
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
      CALL QTRACE(6)
      CALL QUIT('ERROR (GPOPEN) TRYING TO OPEN A NON_EXISTING '//
     &          'FILE AS OLD')
C
C end of GPOPEN
C
      END
C  /* Deck gpclose */
      SUBROUTINE GPCLOSE(LUNIT,DISP)
C
C 23-Feb-2000 K.Ruud
C
C Purpose:
C   General purpose routine for closing files in the Dalton program.
C   The routine will ensure that files that are closed and that will not
C   be reopened in a later part of the program will make their unit 
C   number available for reuse.
C
C   These files are strongly machine dependent, although care has been 
C   taken to avoid using unit numbers that are illegal or reserved on a
C   particular architecture. However, the use of GPOPEN and GPCLSE should
C   remove much of the machine dependence in the rest of the Dalton 
C   program
C
C Input:
C  LUNIT    Mandatory
C  STATUS   Indicate whether the file should be removed or kept
C           ('KEEP' or 'DELETE'). Note that the default is that the file
C           is to be deleted
C
#include <implicit.h>
#include <priunit.h>
      CHARACTER*(*) DISP
      CHARACTER*80 FNNAME, FNNM2
      CHARACTER*3 IODIR
      LOGICAL FEXIST, OPND, OLDDX
      COMMON /UNITAR/ IUNTAB(1:99)
#include <chrnos.h>
C
C     We first deal with the unit number
C
      IF ((LUNIT .LT. 1) .OR. (LUNIT .GT. 99) .OR.
     &    (LUNIT .EQ. 5) .OR. (LUNIT .EQ. 6)) GOTO 9001
ckr     &    (LUNIT .EQ. 5)) GOTO 9001
ckrChj  &    (LUNIT .EQ. 5) .OR. (LUNIT .EQ. 6)) GOTO 9001
ckr      IF (LUNIT .EQ. 6) GOTO 30
C
C     Check that the file actually have been opened
C
      IF (IUNTAB(LUNIT) .EQ. 0) GOTO 9002
C
      INQUIRE(UNIT=LUNIT,EXIST=FEXIST,OPENED=OPND)
      IF (.NOT. FEXIST .OR. .NOT. OPND) THEN
         WRITE (LUPRI,'(/A/A,I3)')
     &     ' WARNING: Tried to close a non-existent or already '//
     &     'closed file', ' Unit number was: ',LUNIT
         CALL FLSHFO(LUPRI)
      END IF
      IF (DISP .EQ. 'KEEP') THEN
         IF (IPRERR.gt.1) write (LUERR,*) 'GPCLOSE keep',LUNIT
         CLOSE (LUNIT,STATUS='KEEP',ERR=9003)
         IUNTAB(LUNIT) = 0
         LUNIT = -10000-LUNIT
         GOTO 30
      ELSE
#if defined (VAR_SPLITFILES)
         INQUIRE(UNIT=LUNIT,DIRECT=IODIR)
         INQUIRE(UNIT=LUNIT,NAME=FNNAME)
         IF (IODIR .EQ. 'YES') THEN
            INQUIRE(UNIT=LUNIT,RECL=LRECL)
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
            LRECL = LRECL/8
#endif
#if !defined (SYS_NEC) && !defined (SYS_IRIX)
            LRECL = LRECL/4
#endif
         END IF
         LN = 1
 15      CONTINUE
         IF (FNNAME(LN:LN) .NE. ' ') THEN
            LN = LN + 1
            GOTO 15
         END IF
         LN = LN - 1
         IF (IODIR .NE. 'YES') THEN
            IF (FNNAME((LN-1):(LN-1)) .EQ. '-') THEN
               CLOSE(LUNIT,STATUS='KEEP')
               LN = LN - 2
               OPEN (UNIT=LUNIT,FILE=FNNAME(1:LN),STATUS='OLD')
            END IF
         END IF
         IF (IPRERR.gt.1) write (LUERR,*) 'GPCLOSE delete',LUNIT
         CLOSE (LUNIT,STATUS='DELETE',ERR=20)
         IUNIT = 0
 16      CONTINUE
         FNNM2 = FNNAME(1:LN)//'-'//CHRNOS(IUNIT)
         LN2 = LN + 2
         INQUIRE(FILE=FNNM2(1:LN2),EXIST=FEXIST)
         IF (FEXIST) THEN
            IF (IODIR .EQ. 'YES') THEN
               CALL OPENDX(LUNIT,FNNM2(1:LN2),LRECL,'OLD',OLDDX)
            ELSE
               OPEN(UNIT=LUNIT,FILE=FNNM2(1:LN2),STATUS='OLD')
            END IF
            CLOSE(LUNIT,STATUS='DELETE')
            IUNIT = IUNIT + 1
            GOTO 16
         END IF
#else
         CLOSE (LUNIT,STATUS='DELETE',ERR=20)
#endif
         IUNTAB(LUNIT) = 0
         LUNIT = -20000-LUNIT
         GOTO 30
      END IF
C
C     We just close it
C     
 20   CLOSE (LUNIT)
      IUNTAB(LUNIT) = 0
      LUNIT = -30000-LUNIT
C
C     We release that saved unit number by resetting it to -(code)-LUNIT
C
 30   CONTINUE
      RETURN
C     
C error branches
C
 9001 CONTINUE
      WRITE (6,'(//A/A,I15/A//A)')
     &   '--> ERROR (GPCLOSE) TRYING TO CLOSE AN ILLEGAL FILE NUMBER',
     &   '--> SOMEBODY HAS TRIED TO CLOSE UNIT',LUNIT,
     &   '--> THE PROGRAM DOES NOT ALLOW THE USE OF THIS UNIT NUMBER',
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
      CALL QTRACE(6)
      CALL QUIT('ERROR (GPCLOSE) ILLEGAL FILE NUMBER REQUESTED')
C
 9002 CONTINUE
      WRITE (6,'(//A/A,I15//A)')
     &   '--> ERROR (GPCLOSE) TRYING TO CLOSE A FILE NOT IN USE',
     &   '--> SOMEBODY IS TRYING TO USE A FILENUMBER THAT HAS '//
     &   'NOT BEEN USED :',LUNIT,
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
      CALL QTRACE(6)
      CALL QUIT('ERROR (GPCLOSE) TRYING TO CLOSE A FILE NOT IN USE')
C
 9003 CONTINUE
      WRITE (6,'(//A/A/A//A)')
     &   '--> ERROR (GPCLOSE) TRYING TO KEEP A SCRATCH FILE',
     &   '--> A FILE HAS BEEN INDICATED TO BE CLOSED AND KEPT,',
     &   '--> BUT IT APPEARS THE FILE IS A SCRATCH FILE',
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
      CALL QTRACE(6)
      CALL QUIT('ERROR (GPCLOSE) TRYING TO CLOSE A SCRATCH FILE ')
C
C end of GPCLOSE
C
      END
C  /* Deck gpinq */
      SUBROUTINE GPINQ(FILENM,TASK,VALUE)
C
C     Generalized routine for checking the existence of a file.
C     Needed for parallel I/O calculations, where the name of the file
C     may have been altered with a processor-specific appendage.
C     The routine is only needed when searching on a filename
C
C     K.Ruud, San Diego Aug 2000
C
#include <implicit.h>
      CHARACTER*(*) FILENM, TASK
      CHARACTER FILTMP*80
      INTEGER FILELN
      LOGICAL VALUE
#include <chrnos.h>
c maxorb for infpar.h
#include <maxorb.h>
#include <infpar.h>
C
      FILELN = LEN(FILENM)
      IF (MYNUM .GT. 0) THEN
         FILTMP(1:FILELN+4) = FILENM(1:FILELN)//
     &                        '#'//CHRNOS(MYNUM/100)//
     &                        CHRNOS(MOD(MYNUM,100)/10)//
     &                        CHRNOS(MOD(MOD(MYNUM,100),10))
         FILELN = FILELN + 4
      ELSE
         FILTMP(1:FILELN) = FILENM(1:FILELN)
      END IF
      IF (TASK(1:5) .EQ. 'EXIST') THEN
         INQUIRE(FILE=FILTMP(1:FILELN),EXIST=VALUE)
      ELSE IF (TASK(1:5) .EQ. 'OPENE') THEN
         INQUIRE(FILE=FILTMP(1:FILELN),OPENED=VALUE)
      ELSE
         CALL QUIT('Unknown inquire status requested')
      END IF
      RETURN
      END
C  /* Deck gpioin */
      SUBROUTINE GPIOIN
C
C 1-Mar-2000 K.Ruud
C 10-Apr-2000 HJAaJ: changed initialization to -9xxx from 0,
C                    because unit 0 is stderr on many systems.
C                    The number -9xxx is used for easier debugging.
C
#include <mxcent.h>
#include <eribuf.h>
#include <priunit.h>
#include <inftap.h>
#include <r12int.h>
C
      LUCMD  = -9001
      LUPRI  = -9002
C
      LUCME  = -9003
      LUMOL  = -9004
      LUPROP = -9005
      LUSOL  = -9006
      LUINTA = -9007
      LUONEL = -9008
      LUSUPM = -9009
      LUTLM  = -9010
      LUDA1  = -9011
      LUITMP = -9012
      LU2DER = -9013
      LUDASP = -9014
      LUAAT  = -9015
      LURDR  = -9016
      LURDI  = -9017
      LUGDR  = -9018
      LUGDI  = -9019
      LUGDT  = -9020
      LURDT  = -9021
      LUDFCK = -9022
      LUSFDA = -9023
      LUFDC  = -9024
      LUWLK  = -9025
      LUPAO  = -9026
      LUPAS  = -9027
      LUNR1  = -9028
      LUNR3  = -9029
      LUNR5  = -9030
      LUCMD  = -9031
      LUINTR = -9032
      LUMOLDEN=-9033
      DO I = 1, MXCOOR
         LUAORC(I) = -9034 - I
      ENDDO
C
C     /R12INT/ (WK/UniKA/26-11-2002).
C
      DO I = 1, 5
         LUR12(I) = -9034 - MXCOOR - I
      ENDDO
C
C     /RSPTAP/
C
      LUAHSO = -9201
      LUCRV1 = -9202
      LUCRV2 = -9203
      LUXYVE = -9204
      LUCRVE = -9205
      LURSP3 = -9206
      LURSP4 = -9207
      LURSP5 = -9208
      LUMHSO = -9209
C
C     /SIRTAP/
C
      LUINTM = -9301
      LUIT1  = -9302
      LUIT2  = -9303
      LUIT3  = -9304
      LUIT5  = -9305
      LUINF  = -9306
      LUH2AC = -9307
      LUSIFC = -9308
C
      FNSOL  = 'AOSOLINT'
      ABARDR = 'ABACUS.RD'
      ABARDI = 'ABACUS.RDI'
      ABAGDR = 'ABACUS.GD'
      ABAGDI = 'ABACUS.GDI'
      ABAGDT = 'ABACUS.GDT'
      ABARDT = 'ABACUS.RDT'
      ABADFK = 'ABACUS.DFK'
      ABASF  = 'ABACUS.SF'
      ABATLM = 'ABACUS.TLM'
      ABAWLK = 'DALTON.WLK'
      ABAIRC = 'DALTON.IRC'
      ABANR1 = 'ABAENR.RST'
      ABANR3 = 'ABAENR.BVC'
      ABANR5 = 'ABAENR.SVC'
      FNINTM = 'MOTWOINT'
      FNSUPM = 'AOSUPINT'
      FNONEL = 'AOONEINT'
      FNSIFC = 'SIRIFC'
      LBSIFC = 'SIR IPH '
C
      RETURN
      END
C /* Deck daintf */
      SUBROUTINE DAINTF(LUNIT,NAME)
C
C     Interface for assigning a unit number to a unit before opening it in
C     DANAME. Needed to ensure compatibility with use of GPOPEN/GPCLOSE
C     K.Ruud, March 28 1999
C
#include <implicit.h>
      CHARACTER NAME*(*)
      COMMON /UNITAR/ IUNTAB(1:99)
C      
      IF (LUNIT .LE. 0) THEN
         IUN = 0
 10      CONTINUE
         IUN = IUN + 1
         IF (IUN .GT. 99) GOTO 9001
         IF ((IUN .EQ. 5) .OR. (IUN .EQ. 6))
     &        GOTO 10
         IF (IUNTAB(IUN) .NE. 0) GOTO 10
         LUNIT = IUN
      END IF
      CALL DANAME(LUNIT,NAME)
      RETURN
C     
C error branch
C
 9001 CONTINUE
      WRITE (6,'(//A/A/A//A)')
     &   '--> ERROR (DAINTF) NO MORE AVAILABLE FILENUMBERS!',
     &   '--> THIS CALCULATION EITHER NEEDS TOO MANY SIMULTANEOUS '//
     &   'FILES OR', 
     &   '--> SOMEBODY HAS FORGOTTEN TO CLOSE FILES IN THE SOURCE CODE',
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
      CALL QTRACE(6)
      CALL QUIT('ERROR (DAINTF) No more available file numbers')
C
      END
