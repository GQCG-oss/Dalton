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
C  /* Deck flshfo */
      SUBROUTINE FLSHFO (IUNIT)
C
C *** THIS SUBROUTINE IS SYSTEM DEPENDENT ***
C
C     Flush formatted output unit (empty buffers).
C     If no flush utility, this is achieved by
C     CLose and reOPen Formatted Output
C
C Written 21-Nov-1983 by Hans Jorgen Aa. Jensen in Uppsala, Sweden.
C Last revision 16-Jul-1984 hjaaj / 30-Oct-1984 hjaaj (extendsize)
C 10-Feb-1989 hjaaj, renamed CLOPFO to FLSHFO
C
C Calls to this subroutine makes it possible to read the output
C up to the moment of the last call while the program continues
C executing (provided the computer allows shared access).
C This subroutine may be a dummy routine.
C
      IF (IUNIT .GE. 0) THEN
C     ... do not try to flush unused units (e.g. LUPRI on a slave) /hjaaj
#if  defined (SYS_AIX)  || defined (SYS_IRIX) || defined (SYS_NEC)  \
  || defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_LINUX) \
  || defined (SYS_SUN)  || defined (SYS_HAL) || defined (SYS_T90)   \
  || defined (SYS_HPUX) || defined (SYS_SX)
C
C        Force transfer of all buffered output to the file or device
C        associated with logical unit IUNIT.
C
         CALL FLUSH(IUNIT)
#endif
      END IF
      RETURN
      END
C  /* Deck getdat */
      SUBROUTINE GETDAT(CDATE,CTIME)
C
C     24-Jan-1988 Hans Joergen Aa. Jensen
C
C     Return date and time as character*8, for labels.
C
      CHARACTER*(8) CDATE, CTIME
#if defined (SYS_CRAY) || defined (SYS_T90) || defined (SYS_T3D)
      INTEGER IDATE, ITIME, DATE, CLOCK
      IDATE = DATE()
      ITIME = CLOCK()
      WRITE (CDATE,'(A8)') IDATE
      WRITE (CTIME,'(A8)') ITIME
#else
#if defined (SYS_AIX) || defined (SYS_DEC) || defined (SYS_IRIX)  \
 || defined (SYS_SUN) || defined (SYS_LINUX) || defined(SYS_HPUX) \
 || defined (SYS_SX) || defined (SYS_DARWIN) || defined (SYS_NEC) 
      CHARACTER*(24) FFDATE
#if defined (SYS_HAL)
      CALL FDATE(FFDATE)
#else
      CHARACTER*(24) FDATE
      FFDATE = FDATE()
#endif
      CDATE  = FFDATE(9:10)//FFDATE(5:7)//FFDATE(23:24)//' '
      CTIME  = FFDATE(12:19)
#else
      CDATE = ' -date- '
      CTIME = ' -time- '
#endif
#endif
      RETURN
      END

      subroutine gettim(cputime,walltime)
      implicit none
      real*8      cputime, walltime

      logical     first
      data        first /.true./

      real*8      TCPU0, TWALL0, tcpu1, twall1
      save        TCPU0, TWALL0
      integer*4   dateandtime0(8), dateandtime1(8)

      if (first) then
         first = .false.
         call cpu_time(TCPU0)
         call date_and_time(values=dateandtime0)
         call get_walltime(dateandtime0,TWALL0)
      end if
      call cpu_time(tcpu1)
      call date_and_time(values=dateandtime1)
      call get_walltime(dateandtime1,twall1)

      cputime  = tcpu1  - TCPU0
      walltime = twall1 - TWALL0

      end subroutine gettim

!> \brief Get elapsed walltime in seconds since 1/1-2010 00:00:00
!> \author S. Host
!> \date October 2010
!>
!> Years that are evenly divisible by 4 are leap years. 
!> Exception: Years that are evenly divisible by 100 are not leap years, 
!> unless they are also evenly divisible by 400. Source: Wikipedia
!>
      subroutine get_walltime(dateandtime,walltime)
      implicit none
      integer*4  dateandtime(8) ! "values" output from fortran intrinsic subroutine date_and_time
      real*8     walltime       ! Elapsed wall time in seconds
      integer*4  month, year
   
! The output from the fortran intrinsic date_and_time
! gives the following values:
! 1. Year
! 2. Month
! 3. Day of the month
! 4. Time difference in minutes from Greenwich Mean Time (GMT)
! 5. Hour
! 6. Minute
! 7. Second
! 8. Millisecond

! Count seconds, minutes, hours, days, months and years and sum up seconds:

      walltime = 1.0d-3*dateandtime(8)                     !Milliseconds
      walltime = walltime + dateandtime(7)*1.0d0           !Seconds counted
      walltime = walltime + 60d0*dateandtime(6)            !Minutes counted
      walltime = walltime + 3600d0*dateandtime(5)          !Hours counted
      walltime = walltime + 24d0*3600d0*(dateandtime(3)-1) !Days counted (substract 1 to count only whole days)

      !Months are special, since they are not equally long:

      do month = 1, dateandtime(2)-1 !substract 1 to count only whole months
         if (month == 1 .or. month == 3 .or. month == 5 .or.
     &       month == 7 .or. month == 8 .or. month == 10) then  !Since we subtract 1, month can never be 12
            walltime = walltime + 31d0*24d0*3600d0
         else if (month == 2) then
            if (.false.) then !insert exception for if current year is a leap year
               walltime = walltime + 29d0*24d0*3600d0
            else
               walltime = walltime + 28d0*24d0*3600d0
            endif
         else if (month == 4 .or. month == 6 .or. month == 9 .or.
     &           month == 11) then
            walltime = walltime + 30d0*24d0*3600d0
         else
            call quit('Unknown month (get_walltime)')
         endif
      enddo

      !Years are special, since leap years are one day longer:

      do year = 2010, dateandtime(1) 
         if (mod(year,400)==0) then
            walltime = walltime + 366*24*3600 !Leap year
         else if (mod(year,100)==0) then
            walltime = walltime + 365*24*3600 !Not leap year
         else if (mod(year,4)==0) then
            walltime = walltime + 366*24*3600 !Leap year
         else
            walltime = walltime + 365*24*3600 !Not leap year
         endif
      enddo

      end subroutine get_walltime
C  /* Deck timtxt */
      SUBROUTINE TIMTXT(TEXT,TIMUSD,LUPRIN)
C
C TIMTXT based on TIMER by TUH //900709-hjaaj
C
#include "implicit.h"
      CHARACTER*(*) TEXT
      CHARACTER AHOUR*6, ASEC*8, AMIN*8
C
      ISECND = NINT(TIMUSD)
      IF (ISECND .GE. 60) THEN
         MINUTE = ISECND/60
         IHOURS = MINUTE/60
         MINUTE = MINUTE - 60*IHOURS
         ISECND = ISECND - 3600*IHOURS - 60*MINUTE
         IF (IHOURS .EQ. 1) THEN
            AHOUR = ' hour '
         ELSE
            AHOUR = ' hours'
         END IF
         IF (MINUTE .EQ. 1) THEN
            AMIN = ' minute '
         ELSE
            AMIN = ' minutes'
         END IF
         IF (ISECND .EQ. 1) THEN
            ASEC = ' second '
         ELSE
            ASEC = ' seconds'
         END IF
         IF (IHOURS .GT. 0) THEN
            WRITE(LUPRIN,100)
     *            TEXT, IHOURS, AHOUR, MINUTE, AMIN, ISECND, ASEC
         ELSE
            WRITE(LUPRIN,200) TEXT, MINUTE, AMIN, ISECND, ASEC
         END IF
      ELSE
         WRITE(LUPRIN,300) TEXT,TIMUSD
      END IF
  100 FORMAT(1X,A,I4,A,I3,A,I3,A)
  200 FORMAT(1X,A,     I3,A,I3,A)
  300 FORMAT(1X,A,F7.2,' seconds')
      RETURN
      END
C  /* Deck tstamp */
      SUBROUTINE TSTAMP(TEXT,LUPRIN)
C
C Copyright Hans Joergen Aa. Jensen 9-Jul-1990
C
C Purpose: To stamp as many as possible of
C          text, date, time, computer, and hostname to LUPRIN
C
#include "implicit.h"
      CHARACTER*(*) TEXT
C
#if defined (SYS_UNIX) || defined (SYS_DEC) || defined (SYS_IRIX) || defined (SYS_NEC) || defined (SYS_SUN) || defined (SYS_LINUX)
      CHARACTER*(40) HSTNAM
      CHARACTER*(24) FDATE
#endif
#if defined (SYS_AIX)
      CHARACTER*(24) fdate
      CHARACTER*(32) HSTNAM
#endif
#if defined (SYS_HAL)
      CHARACTER*(24) CFDATE
      CHARACTER*(32) HSTNAM
#endif
#if defined (SYS_HPUX)
      CHARACTER CFDATE*9, CDATE*8, CTIME*10, HSTNAM*40
#endif
C
      LTEXT = LEN(TEXT)
      IF (LTEXT .GT. 0) THEN
        IF (TEXT(1:LTEXT) .EQ. 'INIT') THEN
CHJ March 2005: this is not working on all machines, so ...
C        WRITE (LUPRIN,'(3(/T3,2A)/)')
C    &   'Last compilation       : ',
C    &   LAST_DALTON_COMPILATION
C    &  ,'Fortran and C compilers: ',
C    &   F_AND_C_COMPILERS
C    &  ,'Scientific libraries   : ',
C    &   LIBRARY_LIST
           WRITE (LUPRIN,'()')
         ELSE
           WRITE (LUPRIN,'(/A)') TEXT
         END IF
      ELSE
         WRITE (LUPRIN,'()')
      END IF
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
      CDATE = DATE()
      CTIME = CLOCK()
      WRITE (LUPRIN,'(T6,A,A8,1X,A8)')
     *   'Date and time (CRAY) : ',CDATE,CTIME
#endif
#if defined (SYS_HAL)
      CALL FDATE(CFDATE)
      WRITE (LUPRIN,'(T6,2A)') 'Date and time (HAL)    : ',CFDATE
#endif
#if defined (SYS_HPUX)
C
C     The following lines have been commented out because they
C     do not work on all HPs
C
C      CALL DATE(CFDATE)
C      CALL DATE_AND_TIME(CDATE,CTIME)
C      WRITE (LUPRIN,'(T6,8A)') 'Date and time (HP-UX)  : ',CFDATE,' ',
C     &     CTIME(1:2),':',CTIME(3:4),':',CTIME(5:6)
C      WRITE (LUPRIN,'(T6,2A)') 'Date and time (HP-UX)  : ',FDATE()
#endif
#if defined (SYS_DEC)
      WRITE (LUPRIN,'(T6,2A)') 'Date and time (DEC)    : ',FDATE()
#endif
#if defined (SYS_IRIX)
      WRITE (LUPRIN,'(T6,2A)') 'Date and time (IRIX)   : ',FDATE()
#endif
#if defined (SYS_NEC)
      WRITE (LUPRIN,'(T6,2A)') 'Date and time (NEC-SX4): ',FDATE()
#endif
#if defined (SYS_SUN)
      WRITE (LUPRIN,'(T6,2A)') 'Date and time (SUN)    : ',FDATE()
#endif
#if defined (SYS_SX)
      WRITE (LUPRIN,'(T6,2A)') 'Data and time (SX)     : ',FDATE()
#endif
#if defined (SYS_LINUX) && !(defined(GFORTRAN) && GFORTRAN <= 402)
c     bug in gfortran-4.0.2 causes a segfault here.
      WRITE (LUPRIN,'(T6,2A)') 'Date and time (Linux)  : ',FDATE()
#endif
#if defined (SYS_LINUX) || defined (SYS_UNIX) || defined (SYS_HPUX)\
 || defined (SYS_DEC) || defined (SYS_IRIX) || defined (SYS_NEC)   \
 || defined (SYS_SUN)
#ifdef VAR_G95
      IRET = HOSTNM(HSTNAM)
#else
      CALL HOSTNM(HSTNAM)
#endif
      WRITE (LUPRIN,'(T6,2A)') 'Host name              : ',HSTNAM
#endif
#if defined (SYS_HAL)
      IRET = HOSTNM(HSTNAM)
      WRITE (LUPRIN,'(T6,2A)') 'Host name              : ',HSTNAM
#endif
#if defined (SYS_AIX)
      WRITE (LUPRIN,'(T6,2A)') 'Date and time (IBM-AIX): ',fdate()
      CALL hostnm(HSTNAM)
      WRITE (LUPRIN,'(T6,2A)') 'Host name              : ',HSTNAM
#endif
C
      RETURN
      END
C  /* Deck ordrss */
      SUBROUTINE ORDRSS(EVEC,EVAL,ISS,N,NEVEC)
C
C 920729-hjaaj (based on ORDER)
C
C Purpose: order the N values in EVAL and their associated vectors
C          in EVEC so EVAL(i+1) .ge. EVAL(i),
C          but only within the class of vectors having the
C          same value in the ISS array (which could be the
C          supersymmetry of the orbital).
C
#include "implicit.h"
      DIMENSION EVEC(*),EVAL(*),ISS(*)
      IF (N.LE.1) RETURN
      IN = 1
      DO 10 I=1,N-1
        EMIN = EVAL(I)
        IMIN = I
        ISSI = ISS(I)
        DO 20 J=I+1,N
C         IF (ISS(J) .NE. ISSI) GO TO 20
C           ... now also reorder diff. supsym, instead update ISS(:)
C               /hjaaj aug 04
          IF (EVAL(J) .LT. EMIN) THEN
            EMIN = EVAL(J)
            IMIN = J
          ENDIF
   20   CONTINUE
        IF (IMIN.NE.I) THEN
          EVAL(IMIN)=EVAL(I)
          EVAL(I)=EMIN
          IF (NEVEC .GT. 0) THEN
            CALL DSWAP(NEVEC,EVEC(IN),1,EVEC((IMIN-1)*NEVEC+1),1)
          ENDIF
          ISS(I) = ISS(IMIN)
          ISS(IMIN) = ISSI
        ENDIF
        IN = IN + NEVEC
   10 CONTINUE
      RETURN
      END
C  /* Deck ord2ss */
      SUBROUTINE ORD2SS(EVEC,EVAL,ISS,N,NEVEC)
C
C Purpose: order the N values in EVAL and their associated vectors
C          in EVEC so EVAL(i+1) .le. EVAL(i) using the infomation in ISS
C          (this is opposite order of "ORDRSS")
C
#include "implicit.h"
      DIMENSION EVEC(*),EVAL(*),ISS(*)
      IF (N.LE.1) RETURN
      IN = 1
      DO 10 I=1,N-1
         EMAX = EVAL(I)
         IMAX = I
         ISSI = ISS(I)
         DO 20 J=I+1,N
C           IF (ISS(J) .NE. ISSI) GO TO 20
C           ... now also reorder diff. supsym, instead update ISS(:)
C               /hjaaj aug 04
            IF (EVAL(J) .GT. EMAX) THEN
               EMAX = EVAL(J)
               IMAX = J
            ENDIF
   20    CONTINUE
         IF (IMAX.NE.I) THEN
            EVAL(IMAX)=EVAL(I)
            EVAL(I)=EMAX
            IF (NEVEC .GT. 0) THEN
              CALL DSWAP(NEVEC,EVEC(IN),1,EVEC((IMAX-1)*NEVEC+1),1)
            ENDIF
            ISS(I) = ISS(IMAX)
            ISS(IMAX) = ISSI
         ENDIF
         IN = IN + NEVEC
   10 CONTINUE
      RETURN
      END
C  /* Deck order */
      SUBROUTINE ORDER(EVEC,EVAL,N,NEVEC)
C
C Purpose: order the N values in EVAL and their associated vectors
C          in EVEC so EVAL(i+1) .ge. EVAL(i)
C
C Revisions:
C   29-Jul-1992 hjaaj (only dswap if nevec .gt. 0)
C    2-Nov-1984 hjaaj (new parameter NEVEC, EVEC(1:NEVEC,1:N))
C   27-Oct-1984 hjaaj (reduced number of swaps)
C
#include "implicit.h"
      DIMENSION EVEC(*),EVAL(*)
      IF (N.LE.1) RETURN
      IN = 1
      DO 10 I=1,N-1
        EMIN = EVAL(I)
        IMIN = I
        DO 20 J=I+1,N
          IF (EVAL(J) .LT. EMIN) THEN
            EMIN = EVAL(J)
            IMIN = J
          ENDIF
   20   CONTINUE
        IF (IMIN.NE.I) THEN
          EVAL(IMIN)=EVAL(I)
          EVAL(I)=EMIN
          IF (NEVEC .GT. 0) THEN
            CALL DSWAP(NEVEC,EVEC(IN),1,EVEC((IMIN-1)*NEVEC+1),1)
          ENDIF
        ENDIF
        IN = IN + NEVEC
   10 CONTINUE
      RETURN
      END
C  /* Deck order2 */
      SUBROUTINE ORDER2(EVEC,EVAL,N,NEVEC)
C
C Purpose: order the N values in EVAL and their associated vectors
C          in EVEC so EVAL(i+1) .le. EVAL(i)
C          (this is opposite order of "ORDER")
C
C Revisions:
C   29-Jul-1992 hjaaj (only dswap if nevec .gt. 0)
C    5-Aug-1985 hjaaj (first version, based on ORDER)
C
#include "implicit.h"
      DIMENSION EVEC(*),EVAL(*)
      IF (N.LE.1) RETURN
      IN = 1
      DO 10 I=1,N-1
         EMAX = EVAL(I)
         IMAX = I
         DO 20 J=I+1,N
            IF (EVAL(J) .GT. EMAX) THEN
               EMAX = EVAL(J)
               IMAX = J
            ENDIF
   20    CONTINUE
         IF (IMAX.NE.I) THEN
            EVAL(IMAX)=EVAL(I)
            EVAL(I)=EMAX
            IF (NEVEC .GT. 0) THEN
              CALL DSWAP(NEVEC,EVEC(IN),1,EVEC((IMAX-1)*NEVEC+1),1)
            ENDIF
         ENDIF
         IN = IN + NEVEC
   10 CONTINUE
      RETURN
      END
C  /* Deck our_own_traceback */
      SUBROUTINE OUR_OWN_TRACEBACK
C
C Written 4-Dec-1983 hjaaj; last revision Jan. 2011
C
#include "implicit.h"
#include "priunit.h"

#if defined (SYS_AIX)
      include 'fexcp.h'
#endif

      INTEGER A,B,C
      SAVE    A,B,C
      DATA    A/1/,B/0/,C/0/

C  920522-hjaaj -- ad hoc routine for creating traceback on IBM-AIX
C  and some other operating systems
C  Note: integer divide by zero is the only error which
C        always will cause an exception
C
 
#if defined (SYS_AIX)
      call SIGNAL(SIGTRAP,xl__trce)
#elif defined (VAR_SECSEC)
      CALL STRACE
#elif defined (SYS_NEC)
      CALL ABORT()
#elif defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
      CALL TRBK()
#else

      C = C + A/B

C The print statements below will onl be printed if the integer divide by zero
C code above does not cause the program to exit.
#endif

Chjaaj apr06:
C     Use stderr on unix/linux systems, if LUPRI not defined yet.
      IF (LUPRI .LT. 0) LUPRI = 0
#if defined (SYS_IRIX)
      WRITE(LUPRI,*) 'IRIX has no obvious system traceback facility.'
#endif
#if defined (SYS_SUN)
      WRITE(LUPRI,*) 'SUN has no obvious system traceback facility.'
#endif
#if defined (SYS_LINUX)
      WRITE(LUPRI,*) 'Linux has no obvious system traceback facility.'
#endif
#if defined (SYS_HPUX)
      WRITE(LUPRI,*) 'HP-UX has no obvious system traceback facility.'
#endif

      RETURN
      END
C  /* Deck canon */
      SUBROUTINE CANON(I,J,K,L)
C     reorder I, J, K, L to canonical 2-electron integral order
#include "implicit.h"
      IP=MAX(I,J)
      JP=I+J-IP
      KP=MAX(K,L)
      LP=K+L-KP
      IF (IP.GT.KP) THEN
         I=IP
         J=JP
         K=KP
         L=LP
      ELSE
         I=KP
         J=LP
         K=IP
         L=JP
         IF(I.NE.K)RETURN
         IF(J.GT.L)RETURN
         J=JP
         L=LP
      END IF
      RETURN
      END
C  /* Deck jaco */
      SUBROUTINE JACO (F,V,NB,NMAX,NROWV,BIG,JBIG)
C
C     F is symmetric packed matrix of dimension NB.
C       The first block of size NMAX will be diagonalized.
C     V is for eigenvectors, only V(NROWV,NMAX) will be referenced.
C       On entry it must correspond the basis vectors corresponding to the
C       F matrix on entry, e.g. unit matrix or AO coefficients for each MO.
C
C Revisions:
C   2-Nov-1984 hjaaj (new parameter NROWV such that
C                     dim(V) = (NROWV,NMAX). This makes
C                     it possible to solve eigenproblem
C                     in a reduced basis but get the
C                     eigenvectors in the original full
C                     basis, e.g. less mo's than ao's)
C  23-Feb-1989 hjaaj  Note that if NROWV = 0 then only
C                     eigenvalues will be calculated,
C                     V matrix will not be referenced.
C  27-Jul-1990 hjaaj  Changed -CX,+SX transformation to +CX,-SX
C                     transformation; probably the -CX,+SX
C                     transformation was responsible for that
C                     the eigenvectors easily changed sign.
C                     Changed initial test on NB. Changed SD.
C                     Optimized IR loop.
C     Jun-1992 ov     Parameters for 0.5, 1.5, ... (for Cray)
C  20-Jul-1992 hjaaj  Changed C1,C2 to THRZER
C  30-oct-1992 hjaaj  zero f(ab) to avoid round-off errors
C                     absolute conv.threshold SD=C1
C  18-aug-2005 wmk    Changed C1 to 1.D-15
#include "implicit.h"
#include "priunit.h"
      DIMENSION F(*),V(NROWV,*)
      DIMENSION BIG(*) ,JBIG(*)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, ROOT2 = 0.707106781186548D0)
      PARAMETER(DP5 = 0.5D0, D1P5 = 1.5D0, D1P375 = 1.375D0,
     *          D3P875 = 3.875D0, DP25 = 0.25D0)
#include "thrzer.h"
      DATA C1,C2,C3,C4,C5,C6/1.D-15,THRZER,1.D-20,1.D-14,1.D-9,1.D-5/
Cwas: DATA C1,C2,C3,C4,C5,C6/THRZER,THRZER,1.D-20,1.D-14,1.D-9,1.D-5/
Cwas: DATA C1,C2,C3,C4,C5,C6/1.D-12,1.D-12,1.D-20,1.D-14,1.D-9,1.D-5/
      IF (NB.LE.1 .OR. NMAX.LE.0) RETURN
Cwas: IF (NB.EQ.1) RETURN !900727-hjaaj
      CALL QENTER('JACO')
      CALL GETTIM(TSTRT, WSTRT)
      DO 190 I=1,NB
         JBIGI=0
         J=MIN(I-1,NMAX)
         IF (J .GT. 0) THEN
            II = (I*I-I)/2
            ABIGI=D0
            DO 18 K=1,J
            IF (ABIGI .GE. ABS(F(II+K))) GO TO  18
               ABIGI=ABS(F(II+K))
               JBIGI=K
   18       CONTINUE
         END IF
         IF (JBIGI .GT. 0) THEN
            JBIG(I) = JBIGI
            BIG(I)  = F(II+JBIGI)
         ELSE
            JBIG(I) = 0
            BIG(I)  = D0
         END IF
  190 CONTINUE
C
#if defined (VAR_OLDCODE)
C 900727-hjaaj:
C SD calculation was done in every Jacobi iteration.
C Now the largest absolute element in F is found once and
C the SD based on that value is used in every iteration.
  410 SD=1.05D 00
      DO 220 J=1,NMAX
         DAB=ABS(F(J*(J+1)/2))
CHJ-861103: commented out next line, it seems to make the loop
C           meaningless (setting SD equal to J=NMAX value always!)
C        IF (SD .GT. DAB) SD=DAB
  220    SD=MAX(SD,DAB)
      SD=MAX(C1,C2*SD)
#else
C 921030-hjaaj: SD = C1 now
      NNB = (NB*NB+NB)/2
C     SD = 1.05D0
C     DO 220 J = 1,NNB
C        SD = MAX(SD, ABS(F(J)) )
C 220 CONTINUE
C     SD=MAX(C1,C2*SD)
      SD=C1
C
      MXITJA = 50*NNB
      ITJACO = 0
  410 ITJACO = ITJACO + 1
      IF (ITJACO .GT. MXITJA) THEN
         CALL QUIT('ERROR: JACO did not converge ...')
      END IF
#endif
      T = D0
      DO 230 I=2,NB
      IF (T .GE.  ABS(BIG(I))) GO TO 230
         T = ABS(BIG(I))
         IB= I
  230 CONTINUE
      IF(T.LT.SD) GO TO 420
         IA =JBIG(IB)
         IAA=IA*(IA-1)/2
         IBB=IB*(IB-1)/2
         DIF=F(IAA+IA)-F(IBB+IB)
         IF( ABS(DIF) .GT. C3) GO TO 271
            SX=ROOT2
            CX=ROOT2
         GO TO 270
  271       T2X2 =BIG(IB)/DIF
            T2X25=T2X2*T2X2
         IF(T2X25 .GT. C4) GO TO 240
            CX=D1
            SX=T2X2
         GO TO 270
  240    IF(T2X25 .GT. C5) GO TO 250
            SX=T2X2*(D1 - D1P5*T2X25)
            CX=D1 - DP5*T2X25
         GO TO 270
  250    IF(T2X25 . GT . C6) GO TO 260
            CX=D1+T2X25*(T2X25*D1P375 - DP5 )
            SX= T2X2*(D1 + T2X25*(T2X25*D3P875 - D1P5))
         GO TO 270
  260       T = DP25  / SQRT(DP25   + T2X25)
            CX= SQRT(DP5   + T)
            SX= SIGN( SQRT(DP5 - T),T2X2)
  270    CONTINUE
#include "ivdep.h"
         DO 275 IR=1,IA
            T        = F(IAA+IR)*SX
            F(IAA+IR)= F(IAA+IR)*CX+F(IBB+IR)*SX
            F(IBB+IR)=-T           +F(IBB+IR)*CX
  275    CONTINUE
         IEAA=IAA+IA
         IEAB=IBB+IA
         TT  =F(IEAB)
         F(IEAB)=BIG(IB)
         IF (JBIG(IA) .EQ. 0) THEN
            IRST = IA   + 1
            IEAR = IEAA + IA
            IEBR = IEAB + 1
         ELSE
            IRST = IA
            IEAR = IEAA
            IEBR = IEAB
         END IF
         DO 390 IR = IRST,NB
#if !defined (VAR_OLDCODE)
            IF (IR .EQ. IA) GO TO 360
C              ... we have checked above that JBIG(IA) .ne. 0
#else
            IF (IR .EQ. IA) THEN
               GO TO 360
C              ... we have checked above that JBIG(IA) .ne. 0
C              IF(JBIG(IR)) 360,380,360
            END IF
#endif
            T      = F(IEAR)*SX
            F(IEAR)= F(IEAR)*CX+F(IEBR)*SX
            F(IEBR)=-T         +F(IEBR)*CX
            T   =F(IEAR)
            IT  =IA
            IF(IR-IB) 340,310,320
  310          F(IEAA)=F(IEAA)*CX+F(IEAB)*SX
C              921030+hjaaj: zero f(ab) to avoid round-off errors
C              F(IEAB)=     TT*CX+F(IEBR)*SX
               F(IEAB)=     D0
               F(IEBR)=    -TT*SX+F(IEBR)*CX
            GO TO 360
  320       IF(ABS(T) .GE.  ABS(F(IEBR))) GO TO 340
               T   =F(IEBR)
               IT  =IB
  340       IF(ABS(T) .LT.  ABS(BIG(IR))) GO TO 350
               BIG(IR)  = T
               JBIG(IR) = IT
            GO TO 380
  350       IF(IA .NE. JBIG(IR) .AND. IB .NE. JBIG(IR))  GO TO 380
  360          K= IEAR - IA
               JBIGI = 0
               IR1=MIN (IR-1,NMAX)
               IF (IR1 .GT. 0) THEN
                  ABIGI = D0
                  DO 370 I=1,IR1
                  IF(ABIGI .GE. ABS(F(K+I)))  GO TO 370
                     ABIGI = ABS(F(K+I))
                     JBIGI =I
  370             CONTINUE
               END IF
               IF (JBIGI .GT. 0) THEN
                  JBIG(IR) = JBIGI
                  BIG(IR)  = F(K+JBIGI)
               ELSE
                  JBIG(IR) = 0
                  BIG(IR)  = D0
               END IF
  380       CONTINUE
               IEAR = IEAR + IR
               IF (IR .GE. IB) THEN
                  IEBR = IEBR + IR
               ELSE
                  IEBR = IEBR + 1
               END IF
  390       CONTINUE
C
#include "ivdep.h"
         DO I=1,NROWV
            VIIA = V(I,IA)
            VIIB = V(I,IB)
            V(I,IA) =  VIIA*CX + VIIB*SX
            V(I,IB) = -VIIA*SX + VIIB*CX
         END DO
      GO TO 410
  420 CONTINUE
c     CALL GETTIM(TEND, WEND)
c     WRITE(LUPRI,'(/A,4I10/A,2F20.2)')
c    &    'JACO -- ITJACO, NB,NMAX,NROWV :',ITJACO, NB,NMAX,NROWV,
c    &    'JACO -- CPUTIME, WALLTIME     :',TEND-TSTRT, WEND-WSTRT
      CALL QEXIT('JACO')
      RETURN
      END
C  /* Deck norm */
      SUBROUTINE NORM(S,VC,N,M,W,THNORM,IRETUR)
C
C revised 14-May-1985 hjaaj (call MPAPV instead of CNTRC)
C revised May 2000 hjaaj (two rounds if small norm)
C
C     COMPUTES SCHMIDT-ORTHONORMALIZED SET OF VECTORS
C         CALLING SEQUENCE PARAMETERS ARE AS FOLLOWS
C              S    METRIC MATRIX STORED AS LOWER TRIANGLE (R*8)
C              VC   LOCATION OF ORIGINAL NON-ORTHONORMAL VECTORS (R*8)
C                   FINAL ORTHONORMALIZED VECTORS REPLACE ORIGINAL SET
C              N    DIMENSION OF BASIS SET (I*4)
C              M    NUMBER OF VECTORS TO BE ORTHONORMALIZED (I*4)
C              W    TEMPORARY WORKING AREA OF 2*N WORDS (R*8)
C         RETURNS
C              NORMAL RETURN ORTHONORMALIZED SET OBTAINED
C              RETURN 1       INITIAL VECTORS AT VC LINEARLY
C                             DEPENDENT WITHIN THRES HOLD (THNORM)
C         AUXILIARY ENTRY
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION S(*), VC(N,M), W(*)
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, THRRND = 0.9D0 )
      IRETUR=0
C
C     N = 1 special case
C
      IF (N .EQ. 1) THEN
         IF (M .EQ. 1) THEN
            IF (VC(1,1)*VC(1,1) .LT. THNORM) THEN
               IRETUR=-1
            ELSE
               VC(1,1) = D1/SQRT(S(1))
            END IF
         END IF
         RETURN
      END IF
C
C     BEGIN OUTERMOST LOOP OVER TRIAL VECTOR SET
C
      DO 20 I=1,M
         ITURN = 0
    1    ITURN = ITURN + 1
         IROUND = 0
C
         CALL MPAPV(N,S,VC(1,I),W)
         TNORM = DDOT(N,VC(1,I),1,W(1),1)
         IF (TNORM .LT. THNORM) THEN
Chj      ... zero vector on input
            IRETUR = -I
            RETURN
         END IF
Chj may2000: normalize input vector
C        (we ignore round-off errors as it is renormalized later)
         TNORM = D1 / SQRT(TNORM)
         CALL DSCAL(N,TNORM,VC(1,I),1)
         CALL MPAPV(N,S,VC(1,I),W)
         TNORM = DDOT(N,VC(1,I),1,W(1),1)
C
C     BEGIN COEFFICIENTS AND NORMALIZATION LOOP
C
         DO 5 J=1,I-1
            T = DDOT(N,VC(1,J),1,W(1),1)
            TNORM = TNORM - T*T
    5       W(N+J) = -T
         IF (TNORM .LT. THNORM) THEN
Chj      ... zero vector after orthogonalization
            IRETUR = I
            RETURN
         END IF
         IF (TNORM .LT. THRRND) IROUND = IROUND + 1
Chj      ... experiments have shown that TNORM as big as
C            0.25 can give a normalization error of 1.0D-7 !
         TNORM = D1/SQRT(TNORM)
         DO J=1,I-1
            W(N+J) = W(N+J)*TNORM
         END DO
         W(N+I) = TNORM
C
C        REPLACE VC(*,I)
C
         CALL DGEMM('N','N',N,1,I,1.D0,
     &              VC,N,
     &              W(N+1),I,0.D0,
     &              W,N)
         CALL DCOPY(N,W,1,VC(1,I),1)
C
         IF (ITURN .EQ. 1 .AND. IROUND .GT. 0) THEN
            IF (IPRSTAT .GT. 0) THEN
               WRITE (LUSTAT,*) 'Info: second round in NORM, I=',I
               CALL QDUMP(LUSTAT)
            END IF
            GO TO 1
         END IF
         IF (IROUND.GT.0) CALL QUIT('NORM: round-off errors, see code')
Chj      ... this ought never to happen ...
   20 CONTINUE
      RETURN
      END
C  /* Deck defopn */
      SUBROUTINE DEFOPN(LU)
C
C 29-Sep-1986 Hans Joergen Aa. Jensen
C
C DEFAULT OPEN :
C
C     OPEN(LU,STATUS='UNKNOWN',FORM='UNFORMATTED') with
C     default file name, following CWBOPN by C.W.Bauschlicher.
C
C     This is especially useful for CRAY machines which do not
C     allow for default names in the OPEN statement (i.e.
C     FILE='name' must be specified if you use OPEN).
C
      REWIND (LU)
      READ   (LU,ERR=1,END=1) A
    1 REWIND (LU)
      RETURN
      END
C  /* Deck readi */
      SUBROUTINE READI (IT,N,INTX)
C
C (30-Jan-1984) hjaaj
C
      DIMENSION INTX(N)
      IF (N .GT. 0) THEN
         READ (IT,END = 10) INTX
      ELSE
         READ (IT,END = 10)
      END IF
      RETURN
   10 CONTINUE
      INTX(N)=-1
      RETURN
      END
C  /* Deck readdi */
      SUBROUTINE READDI(IT,IU,N,IX)
      DIMENSION IX(N)
      READ(IT,REC=IU) IX
      RETURN
      END
C  /* Deck readt */
      SUBROUTINE READT (IT,N,X)
#include "implicit.h"
#include "priunit.h"
      CHARACTER*60 FNAME
      DIMENSION X(N)
      IF (IT .LE. 0) GOTO 30
      READ (IT,END=10,ERR=20) X
      RETURN
 10   CONTINUE
      INQUIRE(UNIT=IT,NAME=FNAME)
      WRITE (LUPRI,*) 'READT: END reading file ',FNAME,
     & ' UNIT ',IT, ' record length ',N
      CALL QUIT('READT: END reading file')
 20   CONTINUE
      INQUIRE(UNIT=IT,NAME=FNAME)
      WRITE (LUPRI,*) 'READT: ERROR reading file ',FNAME,
     & ' UNIT ',IT, ' record length ',N
      CALL QUIT('READT: Error reading file')
 30   CONTINUE
      WRITE (LUPRI,*) 'READT ERRROR: non-positive file unit number ',IT
      CALL QUIT('READT ERRROR: non-positive file unit number')
      END
C  /* Deck writi */
      SUBROUTINE WRITI (IT,N,INTX)
C
C (30-Jan-1984) hjaaj
C
      DIMENSION INTX(N)
      WRITE (IT) INTX
      RETURN
      END
#if defined (VAR_SPLITFILES)
C  /* Deck writsi */
      SUBROUTINE WRITSI(IT,N,INTX,JBUF)
C
C     K.Ruud, April 13 2000
C
#include "implicit.h"
#include "2gbdef.h"
#include "priunit.h"
      DIMENSION INTX(N)
      CHARACTER*80 FNNAME, FNNM2
#include "inftap.h"
#include "chrnos.h"
C
      IF ((JBUF + N) .GT. I2GB) THEN
C
C     Ooops, file will be overfull, need to open a new one......
C
         INQUIRE(UNIT=IT,NAME=FNNAME)
         LN = 1
 10      CONTINUE
         IF (FNNAME(LN:LN) .NE. ' ') THEN
            LN = LN + 1
            GOTO 10
         END IF
         LN = LN - 1
         CALL GPCLOSE(IT,'KEEP')
         I = LN - 1
         IF (FNNAME(I:I) .NE. '-') THEN
            FNNM2 = FNNAME(1:LN)//'-0'
            LN = LN + 2
         ELSE
            READ(FNNAME(LN:),'(I1)') INUM
            INUM = INUM + 1
            IF (INUM .GT. 9) THEN
               WRITE (LUPRI,'(/A)') ' DALTON needs to split a file '//
     &              ' more than 11 times.',
     &              ' This is currently not supported'
               CALL QUIT('Too many splittings of a file')
            END IF
            FNNM2 = FNNAME(1:I)//CHRNOS(INUM)
         END IF
         CALL GPOPEN(IT,FNNM2(1:LN),'UNKNOWN',' ',' ',IDUMMY,
     &               .FALSE.)
         JBUF = 0
      END IF
      WRITE (IT) INTX
      JBUF = JBUF + N
      RETURN
      END
      SUBROUTINE WRITST(IT,N,X,JBUF)
C
C     K.Ruud, April 13 2000
C
#include "implicit.h"
#include "2gbdef.h"
#include "priunit.h"
      DIMENSION X(N)
      CHARACTER*80 FNNAME, FNNM2
#include "inftap.h"
#include "chrnos.h"
C
      WRITE (IT) N
      IF ((JBUF + N + 1) .GT. I2GB) THEN
C
C     Ooops, file will be overfull, need to open a new one......
C
         INQUIRE(UNIT=IT,NAME=FNNAME)
         LN = 1
 10      CONTINUE
         IF (FNNAME(LN:LN) .NE. ' ') THEN
            LN = LN + 1
            GOTO 10
         END IF
         LN = LN - 1
         CALL GPCLOSE(IT,'KEEP')
         I = LN - 1
         IF (FNNAME(I:I) .NE. '-') THEN
            FNNM2 = FNNAME(1:LN)//'-0'
            LN = LN + 2
         ELSE
            READ(FNNAME(LN:),'(I1)') INUM
            INUM = INUM + 1
            IF (INUM .GT. 9) THEN
               WRITE (LUPRI,'(/A)') ' DALTON needs to split a file '//
     &              ' more than 11 times.',
     &              ' This is currently not supported'
               CALL QUIT('Too many splittings of a file')
            END IF
            FNNM2 = FNNAME(1:I)//CHRNOS(INUM)
         END IF
         CALL GPOPEN(IT,FNNM2(1:LN),'UNKNOWN',' ',' ',IDUMMY,
     &               .FALSE.)
         JBUF = 0
      END IF
      WRITE (IT) X
      JBUF = JBUF + N + 1
      RETURN
      END
C  /* Deck readsi */
      SUBROUTINE READSI(IT,N,INTX,JBUF)
C
C     K.Ruud, April 13 2000
C
#include "implicit.h"
#include "2gbdef.h"
#include "priunit.h"
      DIMENSION INTX(N)
      CHARACTER*80 FNNAME, FNNM2
#include "inftap.h"
#include "chrnos.h"
C
      IF ((JBUF + N) .GT. I2GB) THEN
C
C     Ooops, file will be overfull, need to open a new one......
C
         INQUIRE(UNIT=IT,NAME=FNNAME)
         LN = 1
 10      CONTINUE
         IF (FNNAME(LN:LN) .NE. ' ') THEN
            LN = LN + 1
            GOTO 10
         END IF
         LN = LN - 1
         CALL GPCLOSE(IT,'KEEP')
         I = LN - 1
         IF (FNNAME(I:I) .NE. '-') THEN
            FNNM2 = FNNAME(1:LN)//'-0'
            LN = LN + 2
         ELSE
            READ(FNNAME(LN:),'(I1)') INUM
            INUM = INUM + 1
            IF (INUM .GT. 9) THEN
               WRITE (LUPRI,'(/A)') ' DALTON wants to read from a '//
     &              ' file split more than 11 times.',
     &              ' This is currently not supported'
               CALL QUIT('Too many splittings of a file')
            END IF
            FNNM2 = FNNAME(1:I)//CHRNOS(INUM)
         END IF
         CALL GPOPEN(IT,FNNM2(1:LN),'UNKNOWN',' ',' ',IDUMMY,
     &               .FALSE.)
         REWIND (IT)
         JBUF = 0
      END IF
      READ (IT, END = 30) INTX
      JBUF = JBUF + N
      RETURN
 30   INTX(N) = -1
      JBUF = JBUF + 1
      RETURN
      END
      SUBROUTINE READST(IT,N,X,JBUF,READ)
C
C     K.Ruud, Dec 08 2000
C
#include "implicit.h"
#include "2gbdef.h"
#include "priunit.h"
      DIMENSION X(*)
      LOGICAL READ
      CHARACTER*80 FNNAME, FNNM2
#include "inftap.h"
#include "chrnos.h"
C
      READ (IT, END=30) N
      IF ((JBUF + N + 1) .GT. I2GB) THEN
C
C     Ooops, file will be overfull, need to open a new one......
C
         INQUIRE(UNIT=IT,NAME=FNNAME)
         LN = 1
 10      CONTINUE
         IF (FNNAME(LN:LN) .NE. ' ') THEN
            LN = LN + 1
            GOTO 10
         END IF
         LN = LN - 1
         CALL GPCLOSE(IT,'KEEP')
         I = LN - 1
         IF (FNNAME(I:I) .NE. '-') THEN
            FNNM2 = FNNAME(1:LN)//'-0'
            LN = LN + 2
         ELSE
            READ(FNNAME(LN:),'(I1)') INUM
            INUM = INUM + 1
            IF (INUM .GT. 9) THEN
               WRITE (LUPRI,'(/A)') ' DALTON wants to read from a '//
     &              ' file split more than 11 times.',
     &              ' This is currently not supported'
               CALL QUIT('Too many splittings of a file')
            END IF
            FNNM2 = FNNAME(1:I)//CHRNOS(INUM)
         END IF
         CALL GPOPEN(IT,FNNM2(1:LN),'UNKNOWN',' ',' ',IDUMMY,
     &               .FALSE.)
         REWIND (IT)
         JBUF = 0
      END IF
      IF (READ) THEN
         READ (IT, END = 30) (X(I),I = 1, N)
      ELSE
         READ (IT, END = 30)
      END IF
      JBUF = JBUF + N + 1
      RETURN
 30   N = -1
      JBUF = JBUF + 1
      RETURN
      END
#endif
C  /* Deck writdi */
      SUBROUTINE WRITDI(IT,IU,N,IX)
      DIMENSION IX(N)
      WRITE(IT,REC=IU) IX
      RETURN
      END
C  /* Deck writt */
      SUBROUTINE WRITT (IT,N,X)
#include "implicit.h"
      DIMENSION X(N)
      WRITE (IT) X
      RETURN
C
      END
C  /* Deck mollab */
      SUBROUTINE MOLLAB(A,LU,LUERR)
C
C  16-Jun-1986 hjaaj
C  (as SEARCH but CHARACTER*8 instead of REAL*8 labels)
C
C  Purpose:
C     Search for MOLECULE labels on file LU
C
      CHARACTER*8 A, B(4), C
#if defined (VAR_SPLITFILES)
#include "dummy.h"
#endif
      CHARACTER FNNAME*80
      DATA C/'********'/
      IRDERR = 0
#if defined (VAR_SPLITFILES)
      INQUIRE (UNIT=LU,NAME=FNNAME)
      LN = 1
 10   CONTINUE
      IF (FNNAME(LN:LN) .EQ. '-') THEN
         LN = LN - 1
         LUBKP = LU
         CALL GPCLOSE(LU,'KEEP')
         LU = LUBKP
         CALL GPOPEN(LU,FNNAME(1:LN),'OLD','SEQUENTIAL','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         INQUIRE (UNIT=LU,NAME=FNNAME)
         REWIND (LU)
         GOTO 1
      ELSE IF (FNNAME(LN:LN) .EQ. ' ') THEN
         GOTO 1
      END IF
      LN = LN + 1
      GOTO 10
#endif
    1 READ (LU,END=3,ERR=6,IOSTAT=IOSVAL) B
      IRDERR = 0
      IF (B(1).NE.C) GO TO 1
      IF (B(4).NE.A) GO TO 1
      IF (LUERR.LT.0) LUERR = 0
      RETURN
C
    3 IF (LUERR.LT.0) THEN
#if defined (VAR_MFDS)
C 880916-hjaaj -- attempt to fix an IBM problem
C IBM shifts to new file after END= branch (e.g. FTxxF002 instead
C     of FTxxF001), backspace makes LU ready for append.
C 940510-hjaaj: same change for Cray multifile datasets
         BACKSPACE LU
#endif
         LUERR = -1
         RETURN
      ELSE
         INQUIRE (UNIT=LU,NAME=FNNAME)
         WRITE(LUERR,4)A,LU,FNNAME
         CALL QTRACE(LUERR)
         REWIND (LU)
         CALL DMPLAB(LU,LUERR)
         CALL QUIT('ERROR (MOLLAB) MOLECULE label not found on file')
      END IF
    4 FORMAT(/' *** ERROR (MOLLAB), MOLECULE label ',A8,
     *        ' not found on unit',I4/' File name: ',A)
C
    6 IRDERR = IRDERR + 1
      IF (IRDERR .LT. 10) GO TO 1
      IF (LUERR.LT.0) THEN
         LUERR = -2
         RETURN
      ELSE
         WRITE (LUERR,7) LU,A,IOSVAL
         CALL QTRACE(LUERR)
         CALL QUIT('ERROR (MOLLAB) error reading file')
      END IF
    7 FORMAT(/' *** ERROR (MOLLAB), error reading unit',I4,
     *       /T22,'when searching for label ',A8,
     *       /T22,'IOSTAT value :',I7)
      END
C  /* Deck fndlab */
      LOGICAL FUNCTION FNDLAB(A,LU)
C
C 26-May-1985 hjaaj -- logical function version of SEARCH
C 16-Jun-1986 hjaaj -- changed to CHARACTER*8 from REAL*8
C
#include "priunit.h"
      CHARACTER*8 A, B(4), C
      DATA C/'********'/
      IRDERR = 0
    1 READ(LU,END=3,ERR=6)B
      IRDERR = 0
      IF (B(1).NE.C) GO TO 1
      IF (B(4).NE.A) GO TO 1
      FNDLAB = .TRUE.
      GO TO 10
C
    6 IRDERR = IRDERR + 1
      IF (IRDERR .LT. 5) GO TO 1
      GO TO 8
    3 CONTINUE
#if defined (VAR_MFDS)
C 880916-hjaaj -- attempt to fix an IBM problem
C IBM shifts to new file after END= branch (e.g. FTxxF002 instead
C     of FTxxF001), backspace makes LU ready for append.
C 940510-hjaaj: same change for Cray multifile datasets
      BACKSPACE (LU,ERR=7)
C     ... aug07-hjaaj: new ERR branch to avoid problems with backspace on empty files.
      GO TO 8
    7 write (LUPRI,*)
     &  'FNDLAB WARNING: ERR branch for backspace on LU ',LU
      call qdump(LUPRI)
#endif
C
    8 FNDLAB = .FALSE.
C
   10 RETURN
      END
C  /* Deck mollb2 */
      SUBROUTINE MOLLB2(SRCLBL,RTNLBL,LU,LUERR)
C
C  28-Jun-1986 hjaaj
C  (as MOLLAB, but returns two middle labels in RTNLBL(2))
C
C  Purpose:
C     Search for MOLECULE labels on file LU
C
      CHARACTER*8 SRCLBL, RTNLBL(2), B(4), STAR8
      PARAMETER (STAR8 = '********')
C
      IRDERR = 0
    1 READ (LU,END=3,ERR=6,IOSTAT=IOSVAL) B
      IRDERR = 0
      IF (B(1).NE.STAR8)  GO TO 1
      IF (B(4).NE.SRCLBL) GO TO 1
C
         RTNLBL(1) = B(2)
         RTNLBL(2) = B(3)
         IF (LUERR.LT.0) LUERR = 0
         RETURN
C
    3 IF (LUERR.LT.0) THEN
#if defined (VAR_MFDS)
C 880916-hjaaj -- attempt to fix an IBM problem
C IBM shifts to new file after END= branch (e.g. FTxxF002 instead
C     of FTxxF001), backspace makes LU ready for append.
C 940510-hjaaj: same change for Cray multifile datasets
         BACKSPACE LU
#endif
         LUERR = -1
         RETURN
      ELSE
         WRITE (LUERR,4) SRCLBL,LU
         CALL QTRACE(LUERR)
         CALL QUIT('ERROR (MOLLB2) MOLECULE label not found on file')
      END IF
    4 FORMAT(/' *** ERROR (MOLLB2), MOLECULE label ',A8,
     *        ' not found on unit',I4)
C
    6 IRDERR = IRDERR + 1
      IF (IRDERR .LT. 5) GO TO 1
      IF (LUERR.LT.0) THEN
         LUERR = -2
         RETURN
      ELSE
         WRITE (LUERR,7) LU,SRCLBL,IOSVAL
         CALL QTRACE(LUERR)
         CALL QUIT('ERROR (MOLLB2) error reading file')
      END IF
    7 FORMAT(/' *** ERROR (MOLLB2), error reading unit',I4,
     *       /T22,'when searching for label ',A8,
     *       /T22,'IOSTAT value :',I7)
      END
C  /* Deck fndlb2 */
      LOGICAL FUNCTION FNDLB2(SRCLBL,RTNLBL,LU)
C
C  5-Aug-1986 hjaaj
C  (as FNDLAB, but returns two middle labels in RTNLBL(2))
C
      CHARACTER*8 SRCLBL, RTNLBL(2), B(4), STAR8
      PARAMETER (STAR8 = '********')
      IRDERR = 0
    1 READ (LU,END=3,ERR=6) B
      IRDERR = 0
      IF (B(1).NE.STAR8)  GO TO 1
      IF (B(4).NE.SRCLBL) GO TO 1
      FNDLB2    = .TRUE.
      RTNLBL(1) = B(2)
      RTNLBL(2) = B(3)
      GO TO 10
C
    6 IRDERR = IRDERR + 1
      IF (IRDERR .LT. 5) GO TO 1
      GO TO 8
    3 CONTINUE
#if defined (VAR_MFDS)
C 880916-hjaaj -- attempt to fix an IBM problem
C IBM shifts to new file after END= branch (e.g. FTxxF002 instead
C     of FTxxF001), backspace makes LU ready for append.
C 940510-hjaaj: same change for Cray multifile datasets
      BACKSPACE LU
#endif
    8 FNDLB2 = .FALSE.
C
   10 RETURN
      END
C  /* Deck fndlb3 */
      SUBROUTINE FNDLB3(SRCLBL,IVALUE,LU)
C
C  06-Jun-2000 ALig
C  (as FNDLB2, but the IVALUE returns the value of the
C   symmetry representation of the operator SRCLBL and 0 if the
C   label does not exist)
C
#include "implicit.h"
      CHARACTER*8 SRCLBL, B(4), STAR8
      PARAMETER (STAR8 = '********')
      IRDERR = 0
      REWIND(LU)
    1 READ (LU,END=3,ERR=6) B
      IRDERR = 0
      IF (B(1).NE.STAR8)  GO TO 1 
      IF (B(4).NE.SRCLBL) GO TO 1 
      IVALUE    = (ICHAR(B(2)(2:2))-ICHAR('0'))
      GO TO 10
C
    6 IRDERR = IRDERR + 1
      IF (IRDERR .LT. 5) GO TO 1 
      GO TO 8 
    3 CONTINUE
#if defined (VAR_MFDS)
C 880916-hjaaj -- attempt to fix an IBM problem
C IBM shifts to new file after END= branch (e.g. FTxxF002 instead
C     of FTxxF001), backspace makes LU ready for append.
C 940510-hjaaj: same change for Cray multifile datasets
      BACKSPACE LU
#endif
    8 IVALUE = 0
C
   10 CONTINUE
      END   
C
C  /* Deck nxtlab */
      LOGICAL FUNCTION NXTLAB(SRCLBL, RTNLBL, LU)
C
C  3-Nov-1986 hjaaj
C  (find and return next MOLECULE label on LU,
C   NXTLAB false if no label found)
C
      CHARACTER*8 SRCLBL, RTNLBL(2), B(4), STAR8
      PARAMETER ( STAR8 = '********' )
      IRDERR = 0
    1 READ(LU,END=3,ERR=6) B
      IRDERR = 0
      IF (B(1) .NE. STAR8) GO TO 1
      NXTLAB = .TRUE.
      SRCLBL = B(4)
      RTNLBL(1) = B(2)
      RTNLBL(2) = B(3)
      GO TO 10
C
    6 IRDERR = IRDERR + 1
      IF (IRDERR .LT. 5) GO TO 1
      GO TO 8
    3 CONTINUE
#if defined (VAR_MFDS)
C 880916-hjaaj -- attempt to fix an IBM problem
C IBM shifts to new file after END= branch (e.g. FTxxF002 instead
C     of FTxxF001), backspace makes LU ready for append.
C 940510-hjaaj: same change for Cray multifile datasets
      BACKSPACE LU
#endif
    8 NXTLAB = .FALSE.
C
   10 RETURN
      END
C  /* Deck dmplab */
      SUBROUTINE DMPLAB(LU,LUPRI)
C
C 27-Mar-1987 hjaaj -- dump remaining labels on file LU
C
      CHARACTER*8 B(4), C
      PARAMETER ( C = '********' )
C
      WRITE (LUPRI, '(//A,I5)') ' >>> DUMP OF LABELS ON UNIT',LU
      IRDERR = 0
      IREC = 0
    1 READ (LU,END=3,ERR=6,IOSTAT=IOSVAL) B
         IRDERR = 0
         IREC = IREC + 1
         IF (B(1).EQ.C) THEN
            WRITE (LUPRI, '(A,I6,4(2X,A8))') ' Rec. no.',IREC,B
         END IF
      GO TO 1
C
    6 CONTINUE
         IRDERR = IRDERR + 1
         IREC = IREC + 1
         WRITE (LUPRI, '(/A,I6,A,I7)')
     &      ' >>> ERROR reading rec. no.',IREC,'; IOSTAT value',IOSVAL
         IF (IRDERR .LT. 5) GO TO 1
         WRITE (LUPRI, '(/A)')
     &      ' >>> ERROR exit from DMPLAB: 5 consecutive read errors <<<'
    3 CONTINUE
      REWIND (LU)
C
      WRITE (LUPRI, '(/I10,A)') IREC,' records read from file.'
      RETURN
      END
C  /* Deck newlab */
      SUBROUTINE NEWLAB(LABEL,LU,LUERR)
C
C  29-Sep-1988 Hans Joergen Aa. Jensen
C
C  Write new MOLECULE-type label to LU
C
      CHARACTER*8 LABEL, LTIME, LDATE, LSTARS
      DATA LSTARS /'********'/
      CALL GETDAT(LDATE,LTIME)
      WRITE (LU,ERR=1000,IOSTAT=IOSVAL) LSTARS,LDATE,LTIME,LABEL
      RETURN
C
 1000 IF (LUERR .LT. 0) THEN
         LUERR = -2
         RETURN
      ELSE
         WRITE (LUERR,'(/3A,I5/A,I7)')
     &      ' NEWLAB: error writing label "',LABEL,'" to unit',LU,
     &      '         IOSTAT value',IOSVAL
         CALL QTRACE(LUERR)
         CALL QUIT('NEWLAB: output error writing label')
      END IF
      END
C  /* Deck newlb2 */
      SUBROUTINE NEWLB2(LABEL,RTNLBL,LU,LUERR)
C
C  29-Sep-1988 Hans Joergen Aa. Jensen
C
C  Write new MOLECULE-type label to LU
C
      CHARACTER*8 LABEL, RTNLBL(2), LSTARS
      DATA LSTARS /'********'/
      WRITE (LU,ERR=1000,IOSTAT=IOSVAL) LSTARS,RTNLBL,LABEL
      RETURN
C
 1000 IF (LUERR .LT. 0) THEN
         LUERR = -2
         RETURN
      ELSE
         WRITE (LUERR,'(/3A,I5/A,I7)')
     &      ' NEWLB2: error writing label "',LABEL,'" to unit',LU,
     &      '         IOSTAT value',IOSVAL
         CALL QTRACE(LUERR)
         CALL QUIT('NEWLB2: output error writing label')
      END IF
      END
C  /* Deck second */
#if !defined (SYS_CRAY) && !defined (SYS_T90) && !defined (VAR_GFORTRAN)
      REAL*8 FUNCTION SECOND ()
      REAL*4 TIME_CPU
      CALL CPU_TIME(TIME_CPU)
      SECOND = TIME_CPU
      RETURN
      END
#endif
C  /* Deck sotmat */
      SUBROUTINE SOTMAT(NMO,UMO,IFAIL)
C
C  16-Feb-1986 Hans Jorgen Aa. Jensen
C
C  Purpose:
C
C    Construct the orbital transformation matrix
C    for transforming a CAS CI vector using a sequence
C    of single orbital transformations as described by
C    Per-Ake Malmquist.
C
C    The matrix is
C
C      C   +  C   = (1 - L) + U(inv);
C       L      U
C
C    where L and U constitute the LU decomposition of the
C    orthogonal orbital transformation matrix UMO.
C
#include "implicit.h"
      DIMENSION UMO(NMO,NMO)
#include "priunit.h"
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0 )
      PARAMETER ( THRESH = 1.D-4 )
C
      CALL QENTER('SOTMAT')
C
C     STEP 1: The LU decomposition
C     ============================
C
      DO 220 K = 1,NMO
         X = UMO(K,K)
         IF (ABS(X) .LE. THRESH) GO TO 960
         X = D1 / X
         DO 120 I = K+1,NMO
            UMO(I,K) = UMO(I,K) * X
  120    CONTINUE
         DO 200 J = (K+1),NMO
            Y = UMO(K,J)
            DO 180 I = K+1,NMO
               UMO(I,J) = UMO(I,J) - UMO(I,K)*Y
  180       CONTINUE
  200    CONTINUE
  220 CONTINUE
C
C     STEP 2: U inverse
C     =================
C
      DET = D1
      DO 300 K = 1,NMO
         DET = DET * UMO(K,K)
  300 CONTINUE
#if defined (VAR_SOTMATTEST)
C
C
      WRITE (LUPRI,'(//A,1P,D10.2)') ' SOTMAT: UMO determinant =',DET
      WRITE (LUPRI,'(//A)') ' SOTMAT: LU matrix'
      CALL OUTPUT (UMO,1,NMO,1,NMO,NMO,NMO,1,6)
#endif
      IF (ABS(DET) .LE. THRESH) GO TO 980
C
      DO 400 K = 1,NMO
         UIDIAG = D1 / UMO(K,K)
         UMO(K,K) = UIDIAG
         DO 380 J = 1,(K-1)
            SUM = D0
            DO 360 I = J,(K-1)
               SUM = SUM - UMO(J,I)*UMO(I,K)
  360       CONTINUE
            UMO(J,K) = UIDIAG * SUM
  380    CONTINUE
  400 CONTINUE
C
C     STEP 3: construct 1 - L
C     =======================
C
      DO 600 K = 1,NMO
         DO 580 J = (K+1),NMO
            UMO(J,K) = -UMO(J,K)
  580    CONTINUE
  600 CONTINUE
C
C     NORMAL AND ERROR RETURNS
C     ========================
C
      IFAIL = 0
      GO TO 9999
C
  960 CONTINUE
      WRITE (LUPRI,'(///A,1P,D10.2/)')
     1   ' >>> SOTMAT: DIAGONAL ELEMENT TOO SMALL:',X
      IFAIL = 1
      GO TO 9999
C
  980 CONTINUE
      WRITE (LUPRI,'(///2A,1P,D10.2/)')
     1   ' >>> SOTMAT: U MATRIX TOO CLOSE TO SINGULARITY,',
     2   ' DETERMINANT =',DET
      IFAIL = 2
C
 9999 CALL QEXIT('SOTMAT')
      RETURN
C
C     END OF SOTMAT.
C
      END
C  /* Deck nofdia */
      INTEGER FUNCTION NOFDIA(N,NDIM,AMAT,THREQL)
C
C  2-Jul-1992 hjaaj
C  This function returns the number of off-diagonal elements
C  with absolute value greater than THREQL.
C
#include "implicit.h"
      DIMENSION AMAT(NDIM,NDIM)
C
      CALL QENTER('NOFDIA')
      NELEM = 0
      DO 220 K = 1,N
         DO 210 I = 1,N
            IF (I .NE. K .AND. ABS(AMAT(I,K)) .GT. THREQL) THEN
               NELEM = NELEM + 1
            END IF
  210    CONTINUE
  220 CONTINUE
      NOFDIA = NELEM
      CALL QEXIT('NOFDIA')
      RETURN
      END
#if defined (SYS_HPUX)
      character*24 function fdate()
      fdate = '** no date available **'
      return
      end
#endif
C  /* Deck matsym */
      INTEGER FUNCTION MATSYM(N,NDIM,AMAT,THRZER)
C
C  27-Apr-1999 hjaaj
C  This function returns 
C   1 if the matrix is symmetric to threshold THRZER
C   2 if the matrix is antisymmetric to threshold THRZER
C   3 if all elements are below THRZER
C   0 otherwise (the matrix is unsymmetric about the diagonal)
C
#include "implicit.h"
      DIMENSION AMAT(NDIM,NDIM)
C
      CALL QENTER('MATSYM')
      ISYM  = 1
      IASYM = 2
      DO 220 J = 1,N
         DO 210 I = 1,J
            AMATS = ABS(AMAT(I,J) + AMAT(J,I))
            AMATA = ABS(AMAT(I,J) - AMAT(J,I))
            IF (AMATS .GT. THRZER) IASYM = 0
C           ... i.e. not antisymmetric
            IF (AMATA .GT. THRZER) ISYM = 0
C           ... i.e. not symmetric
  210    CONTINUE
  220 CONTINUE
      MATSYM = ISYM + IASYM
      CALL QEXIT('MATSYM')
      RETURN
      END
C  /* Deck rewmot */
      SUBROUTINE REWSPL(LU)
C
C     Short interface routine for rewinding a file that may have been 
C     split, to ensure that we search for labels on the first of the split
C     files. This routines preserves the UNIT number.
C
C     K.Ruud, April 19 (2000)
C
#if !defined (VAR_SPLITFILES)
#include "priunit.h"
      INTEGER LU
      IF (LU .LT. 1) THEN
         WRITE (LUPRI,'(/A,I10)')
     &      'ERROR in REWSPL: negative unit number =',LU
         CALL QUIT('REWSPL called with negative unit number LU')
      END IF
      REWIND (LU)
#else
#include "implicit.h"
#include "dummy.h"
      CHARACTER FNNAME*80
C
      INQUIRE (UNIT=LU,NAME=FNNAME)
      LN = 1
 14   CONTINUE
      IF (FNNAME(LN:LN) .EQ. '-') THEN
         LN = LN - 1
         LUBKP = LU
         CALL GPCLOSE(LU,'KEEP')
         LU = LUBKP
         CALL GPOPEN(LU,FNNAME(1:LN),'OLD','SEQUENTIAL',
     &               'UNFORMATTED',IDUMMY,.FALSE.)
         INQUIRE (UNIT=LU,NAME=FNNAME)
         GOTO 15
      ELSE IF (FNNAME(LN:LN) .EQ. ' ') THEN
         GOTO 15
      END IF
      LN = LN + 1
      GOTO 14
 15   CONTINUE
      REWIND (LU)
#endif
      RETURN
      END
C  /* Deck fndmin */
      SUBROUTINE FNDMIN(NELMNT,IPLACE,VEC,NVEC,WRK,LWRK)
C
C 23-Nov-2000 hjaaj (FNDMIN = CIFNMN from 12-Aug-1990 hjaaj)
C (Find minimum elemnts)
C
C Purpose:
C   Return in IPLACE addresses on lowest NELMNT elements in VEC.
C
#include "implicit.h"
      DIMENSION VEC(NVEC), IPLACE(NELMNT), WRK(LWRK)
C
      IF (LWRK .LT. NELMNT) THEN
         CALL QUIT('FNDMIN: Insufficient memory (LWRK .lt. NELMNT)')
      END IF
      IF (NELMNT .GT. NVEC) THEN
         CALL QUIT('FNDMIN ERROR: NELMNT .gt. NVEC')
      END IF
C
C     Sort first NELMNT elements of VEC
C
      DO 120 I = 1,NELMNT
         VECI = VEC(I)
         DO 115 J = 1,(I-1)
         IF (WRK(J) .GT. VECI) THEN
            DO 112 K = I,(J+1),-1
               WRK(K)    = WRK(K-1)
               IPLACE(K) = IPLACE(K-1)
  112       CONTINUE
            IPLACE(J) = I
            WRK(J)    = VECI
         GO TO 120
         END IF
  115    CONTINUE
         IPLACE(I) = I
         WRK(I)    = VECI
  120 CONTINUE
C
C     Find lowest elements by insertion sort
C
      XHGH = WRK(NELMNT)
      DO 140 I = NELMNT+1,NVEC
      IF (VEC(I).GE.XHGH) GO TO 140
         DO 130 J = 1,NELMNT
         IF (VEC(I) .LT. WRK(J)) THEN
            DO 135 K = NELMNT,(J+1),-1
               WRK(K) = WRK(K-1)
               IPLACE(K) = IPLACE(K-1)
  135       CONTINUE
            IPLACE(J) = I
            WRK(J) = VEC(I)
            XHGH   = WRK(NELMNT)
            GO TO 140
         END IF
  130    CONTINUE
  140 CONTINUE
      RETURN
      END
      SUBROUTINE COMPOFF(WRK, WORK, KBASE)
C
C     Nov. 2004 Hans Joergen Aa. Jensen; revised Dec. 2010
C
C     COMPOFF - Compute off-set from WORK to WRK
C               i.e. WRK(I) = WORK(KBASE+I)
C
      REAL*8                    :: WRK(*), WORK(*)
      integer(8), intent(inout) :: kbase
      integer(8)                :: kaddr

      !print *,'LOC(WRK) ',LOC(WRK)
      !print *,'LOC(WORK)',LOC(WORK)
      KADDR = LOC(WRK)  - LOC(WORK)
      KBASE = KADDR / 8 + 1 ! from byte address to REAL*8 address

      !print *,'COMPOFF: KADDR, KBASE', KADDR, KBASE, KBASE*8
      END
C -- end of gphjj.F --
