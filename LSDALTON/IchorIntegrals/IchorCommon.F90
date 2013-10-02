!> @file
!> Contains common routine for the Ichor Code
!> \brief Contains common routine for the Ichor Code
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorCommonModule
use IchorprecisionModule
public:: IchorQuit, ichor_tstamp, Ichor_gettim, ichor_get_walltime, &
     & ichor_timtxt, ichortimer
private
CONTAINS
subroutine ichor_tstamp(TEXT,LUPRIN)
  implicit none
  CHARACTER(*), intent(in) :: TEXT
  integer, intent(in)      :: luprin
  integer :: ltext
#if  defined (SYS_LINUX)
  CHARACTER(40) :: HSTNAM
  CHARACTER(24) :: FDATE
#endif
#if defined (SYS_AIX)
  CHARACTER(24) :: fdate
  CHARACTER(32) :: HSTNAM
#endif    
  LTEXT = LEN(TEXT)
  IF (LTEXT .GT. 0) THEN
     WRITE (LUPRIN,'(/A)') TEXT
  ELSE
     WRITE (LUPRIN,'()')
  END IF
#if defined (SYS_LINUX)
  WRITE (LUPRIN,'(T6,2A)') 'Date and time (Linux)  : ',FDATE()
  CALL HOSTNM(HSTNAM)
  WRITE (LUPRIN,'(T6,2A)') 'Host name              : ',HSTNAM
#endif
#if defined (SYS_AIX)
  ! 930414-hjaaj: apparent error IBM's xlf library routines:
  !    when 'T6' then column 1-5 not blanked but contains text
  !    from a previous print statement!
  !     AIX XL FORTRAN version 2.3+
  WRITE (LUPRIN,'(2A)') '     Date and time (IBM-AIX): ',fdate()
  !     WRITE (LUPRIN,'(T6,2A)') 'Date and time (IBM-AIX): ',fdate_()
  !     CALL hostnm_(HSTNAM)
  CALL hostnm(HSTNAM)
  WRITE (LUPRIN,'(2A)') '     Host name              : ',HSTNAM
  !     WRITE (LUPRIN,'(T6,2A)') 'Host name              : ',HSTNAM
#endif
end subroutine ichor_tstamp

!Return elapsed CPU time and elapsed real time.
subroutine Ichor_gettim(cputime,walltime)
  implicit none
  real(realk), intent(out) :: cputime, walltime  
  real(realk),PARAMETER :: D0 = 0.0E0_realk
  logical    :: first = .true.
  real(realk), save :: TCPU0, twall0
  real(realk)       :: tcpu1, twall1
  integer           :: dateandtime0(8), dateandtime1(8)  
  if (first) then
     first = .false.
     call cpu_time(TCPU0)
     call date_and_time(values=dateandtime0)
     call ichor_get_walltime(dateandtime0,twall0)
  end if
  call cpu_time(tcpu1)
  call date_and_time(values=dateandtime1)
  call ichor_get_walltime(dateandtime1,twall1)
  cputime = tcpu1 - TCPU0
  walltime = twall1 - twall0
end subroutine Ichor_gettim

!> \brief Get elapsed walltime in seconds since 1/1-2010 00:00:00
!> \author S. Host
!> \date October 2010
!>
!> Years that are evenly divisible by 4 are leap years. 
!> Exception: Years that are evenly divisible by 100 are not leap years, 
!> unless they are also evenly divisible by 400. Source: Wikipedia
!>
subroutine ichor_get_walltime(dateandtime,walltime)
implicit none
!> "values" output from fortran intrinsic subroutine date_and_time
integer, intent(in) :: dateandtime(8)
!> Elapsed wall time in seconds
real(realk), intent(out) :: walltime
integer                  :: month, year

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
! We don't count milliseconds.
walltime = 0.0E0_realk
walltime = dble(dateandtime(8))/1.0E3_realk                !Seconds counted
walltime = walltime + dble(dateandtime(7))                 !Seconds counted
walltime = walltime + 60E0_realk*dateandtime(6)            !Minutes counted
walltime = walltime + 3600E0_realk*dateandtime(5)          !Hours counted
walltime = walltime + 24E0_realk*3600E0_realk*(dateandtime(3)-1) !Days counted (substract 1 to count only whole days)

!Months are special, since they are not equally long:
do month = 1, dateandtime(2)-1 !substract 1 to count only whole months
 if (month == 1 .or. month == 3 .or. month == 5 .or. month == 7 .or. &
      &  month == 8 .or. month == 10) then !Since we subtract 1, month can never be 12
    walltime = walltime + 31E0_realk*24E0_realk*3600E0_realk
 else if (month == 2) then
  if (.false.) then !insert exception for if current year is a leap year
     walltime = walltime + 29E0_realk*24E0_realk*3600E0_realk
  else
     walltime = walltime + 28E0_realk*24E0_realk*3600E0_realk
  endif
 else if (month == 4 .or. month == 6 .or. month == 9 .or. month == 11) then
    walltime = walltime + 30E0_realk*24E0_realk*3600E0_realk
 else
    stop 'Unknown month (get_walltime)'
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
end subroutine ichor_get_walltime

subroutine ichor_timtxt(TEXT,TIMUSD,LUPRIN)
  implicit none 
  CHARACTER(*),intent(in) :: TEXT
  real(realk), intent(in) :: timusd
  integer, intent(in)     :: luprin
  CHARACTER :: AHOUR*6, ASEC*8, AMIN*8
  integer :: minute, isecnd, ihours  
  ISECND = NINT(TIMUSD)
  IF (ISECND > 60) THEN
     MINUTE = ISECND/60
     IHOURS = MINUTE/60
     MINUTE = MINUTE - 60*IHOURS
     ISECND = ISECND - 3600*IHOURS - 60*MINUTE
     IF (IHOURS == 1) THEN
        AHOUR = ' hour '
     ELSE
        AHOUR = ' hours'
     END IF
     IF (MINUTE == 1) THEN
        AMIN = ' minute '
     ELSE
        AMIN = ' minutes'
     END IF
     IF (ISECND == 1) THEN
        ASEC = ' second '
     ELSE
        ASEC = ' seconds'
     END IF
     IF (IHOURS /= 0) THEN
        WRITE(LUPRIN,"(1X,A,I4,A,I3,A,I3,A)") TEXT, IHOURS, AHOUR, MINUTE, AMIN, ISECND, ASEC
     ELSE
        WRITE(LUPRIN,"(1X,A,     I3,A,I3,A)") TEXT, MINUTE, AMIN, ISECND, ASEC
     END IF
  ELSE
     WRITE(LUPRIN,"(1X,A,F7.2,' seconds')") TEXT,TIMUSD
  END IF
  
  CALL FLUSH(LUPRIN)
end subroutine ichor_timtxt

subroutine IchorQuit(text,lupri)
implicit none
!> Text string to be printed
CHARACTER(len=*), intent(in) :: TEXT
!> Logical unit number for output
integer, intent(in) :: lupri
integer             :: luprin
real(realk) :: CTOT,WTOT
!
!     Stamp date and time and hostname to output
!
if (lupri > 0) then
   luprin = lupri
else
   luprin = 6
endif

CALL ichor_TSTAMP('  --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---', LUPRIN)
WRITE (LUPRIN,'(/2A/)') ' Reason: ',TEXT

#if  defined (SYS_AIX) || defined (SYS_LINUX)
!     Write to stderr
WRITE (0,'(/A/1X,A)') '  --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---',TEXT
#endif

CALL Ichor_GETTIM(CTOT,WTOT)
CALL Ichor_TIMTXT('>>>> Total CPU  time used in DALTON:',CTOT,LUPRIN)
CALL Ichor_TIMTXT('>>>> Total wall time used in DALTON:',WTOT,LUPRIN)
CALL FLUSH(LUPRIN)
!TRACEBACK INFO TO SEE WHERE IT CRASHED!!
#ifdef VAR_IFORT
CALL TRACEBACKQQ()
#endif
#if defined (SYS_LINUX)
CALL EXIT(100)
#else
STOP 100
#endif
end subroutine IchorQuit

!> \brief take time
!> \author T. Kjaergaard
!> \date 2010
!> \param text label to print along with timings
!> \param CPUTIME the cpu time
!> \param WALLTIME the wall time
!> \param lupri the logical unit number 
SUBROUTINE IchorTimer(TEXT,CPUTIME,WALLTIME,LUPRI,FORCEPRINT)
implicit none
logical,optional  :: FORCEPRINT
INTEGER           :: LUPRI,length
CHARACTER*(*)     :: TEXT
CHARACTER(len=30) :: PRINTTEXT
REAL(REALK)       :: TIME1,TIME2,DELTAWALL,CPUTIME,WALLTIME,DELTACPU

INTEGER :: IUNIT
logical :: FORCEPRINT2
FORCEPRINT2 = .FALSE.
IF(present(FORCEPRINT))THEN
   FORCEPRINT2 = FORCEPRINT
ENDIF

IF (LUPRI.EQ.-1) THEN
   IUNIT = 6
ELSE
   IUNIT = LUPRI
ENDIF

length = LEN(TEXT)
IF(length .GT. 30)THEN
   CALL IchorQuit('TEXTLENGTH PROVIDED TO ICHORTIMER IS LIMITED TO 30',lupri)
ENDIF 
PRINTTEXT='               '
PRINTTEXT(1:length) =  TEXT(1:length) 
IF (PRINTTEXT(1:5) .EQ. 'START') THEN
   CALL ICHOR_GETTIM(CPUTIME,WALLTIME)
ELSE
   CALL ICHOR_GETTIM(TIME1,TIME2)
   DELTACPU=TIME1-CPUTIME
   DELTAWALL=TIME2-WALLTIME
   IF(FORCEPRINT2) THEN
      CALL ICHOR_TIMTXT('>>>  CPU Time used in '//PRINTTEXT//' is',DELTACPU,IUNIT)
      CALL ICHOR_TIMTXT('>>> wall Time used in '//PRINTTEXT//' is',DELTAWALL,IUNIT)
   ENDIF
   WALLTIME=TIME2
   CPUTIME=TIME1
   CALL FLUSH(IUNIT)
ENDIF

END SUBROUTINE ICHORTIMER

END MODULE IchorCommonModule
