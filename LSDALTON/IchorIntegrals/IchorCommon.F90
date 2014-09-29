!> @file
!> Contains common routine for the Ichor Code
!> \brief Contains common routine for the Ichor Code
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorCommonModule
use IchorprecisionModule
!public:: IchorQuit, ichor_tstamp, Ichor_gettim, ichor_get_walltime, &
!     & ichor_timtxt, ichortimer, GenerateOrderdListOfTypes
!private
public
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
#ifdef VAR_IFORT
  use IFCORE
#endif
implicit none
!> Text string to be printed
CHARACTER(len=*), intent(in) :: TEXT
!> Logical unit number for output
integer, intent(in) :: lupri
integer             :: luprin,qqstatus,user_exit_code
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
CALL Ichor_TIMTXT('>>>> Total CPU  time used in LSDALTON:',CTOT,LUPRIN)
CALL Ichor_TIMTXT('>>>> Total wall time used in LSDALTON:',WTOT,LUPRIN)
CALL FLUSH(LUPRIN)
#ifdef VAR_IFORT
#ifndef VAR_INT64
!TRACEBACK INFO TO SEE WHERE IT CRASHED!!
qqstatus=-1
user_exit_code = -1
CALL TRACEBACKQQ("Ichor Called TraceBack:",user_exit_code,qqstatus)
#endif
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

subroutine GenerateOrderdListOfTypes(lupri,nTypes,AngmomOfType,OrderdList)
implicit none
!> nTypes is the number of different types of shells, each type is defined by 
!> an angular momentum, a number of primitives(nPrim), a number of contracted functions
!> (nCont) a set of exponents and a set of contraction coefficients. 
integer,intent(in) :: nTypes
!> AngmomOfType is the angular momentum for each type. 
Integer,intent(in) :: AngmomOfType(ntypes)
!> Orderd List
Integer :: OrderdList(ntypes)
!> Logical unit number of output file.
Integer :: lupri
!local variables
integer :: ILIST,Angmom,Itype

ILIST = 1
!I want the highest angular momentum first
DO Angmom=10,0,-1
  DO Itype=1,nTypes
   IF(Angmom .EQ. AngmomOfType(Itype))THEN
      OrderdList(ILIST) = Itype
      ILIST = ILIST + 1
   ENDIF
  ENDDO
ENDDO
IF(ILIST-1.NE.nTypes)call ichorQuit('GenerateOrderdListOfTypes Error',-1)
end subroutine GenerateOrderdListOfTypes

subroutine build_exp_ContractCoeff_center(nPrimA,nContA,nAtomsA,ntypesA,iTypeA,&
     & exponentsOfTypeA,ContractCoeffOfTypeA,Acenters,startOrbitalOfTypeA,&
     & expA,ContractCoeffA,Acenter,startOrbitalA,MaxnAtomsA,MaxnprimA,MaxnContA,lupri)
  implicit none
  integer,intent(in) :: nPrimA,nContA,nAtomsA,ntypesA,iTypeA,MaxnPrimA,MaxnContA,lupri,MaxnAtomsA
  real(realk),intent(in) :: exponentsOfTypeA(MaxnprimA,ntypesA)
  real(realk),intent(in) :: ContractCoeffOfTypeA(MaxnprimA,MaxnContA,ntypesA)
  real(realk),intent(in) :: Acenters(3,MaxnAtomsA,ntypesA)
  integer,intent(in) :: startOrbitalOfTypeA(MaxnAtomsA,ntypesA)
  real(realk),intent(inout) :: expA(nPrimA),ContractCoeffA(nPrimA,nContA),Acenter(3,nAtomsA)
  integer,intent(inout) :: startOrbitalA(nAtomsA)
  !local variables
  integer :: iPrimA,iContA,iAtomA
  do iPrimA = 1, nPrimA
     expA(iPrimA) = exponentsOfTypeA(iPrimA,iTypeA)
  enddo
  do iContA = 1, nContA
     do iPrimA = 1, nPrimA
        ContractCoeffA(iPrimA,iContA) = ContractCoeffOfTypeA(iPrimA,iContA,iTypeA)
     enddo
  enddo
  do iAtomA = 1, nAtomsA
     Acenter(1,iAtomA) = Acenters(1,iAtomA,itypeA)
     Acenter(2,iAtomA) = Acenters(2,iAtomA,itypeA)
     Acenter(3,iAtomA) = Acenters(3,iAtomA,itypeA)
     startOrbitalA(iAtomA) = startOrbitalOfTypeA(iAtomA,itypeA)
  enddo
end subroutine build_exp_ContractCoeff_center

subroutine build_expP(nPrimA,nPrimB,expA,expB,expP)
  implicit none
  integer,intent(in) :: nPrimA,nPrimB
  real(realk),intent(in) :: expA(nPrimA),expB(nPrimB)
  real(realk),intent(inout) :: expP(nPrimA,nPrimB)
  !local variables
  integer :: i1,i2
  real(realk) :: e1,e2 
  DO i2=1,nPrimB
     e2  = expB(i2)       
     DO i1=1,nPrimA
        e1  = expA(i1)
        expP(i1,i2) = e1 + e2
     ENDDO
  ENDDO
end subroutine build_expP

subroutine build_expQ_inverseexpQ(nPrimC,nPrimD,expC,expD,expQ,inversexpQ)
  implicit none
  integer,intent(in) :: nPrimC,nPrimD
  real(realk),intent(in) :: expC(nPrimC),expD(nPrimD)
  real(realk),intent(inout) :: expQ(nPrimC,nPrimD),inversexpQ(nPrimC,nPrimD)
  !local variables
  integer :: i1,i2
  real(realk) :: e1,e2
  DO i2=1,nPrimD
     e2  = expD(i2)       
     DO i1=1,nPrimC
        e1  = expC(i1)
        expQ(i1,i2) = e1 + e2
        inversexpQ(i1,i2) = 1.0E0_realk/(e1 + e2)
     ENDDO
  ENDDO
end subroutine build_expQ_inverseexpQ

!!$subroutine build_reducedExponents_integralPrefactorPQ(nPrimP,nPrimQ,expQ,expP,&
!!$     & reducedExponents,integralPrefactor)
!!$  implicit none      
!!$  integer,intent(in) :: nPrimP,nPrimQ
!!$  real(realk),intent(in) :: expQ(nPrimQ),expP(nPrimP)
!!$  real(realk),intent(inout) :: reducedExponents(nPrimP,nPrimQ)
!!$  real(realk),intent(inout) :: integralPrefactor(nPrimP,nPrimQ)
!!$  !local variables 
!!$  integer :: iPrimQ,iPrimP
!!$  real(realk) :: p,q,p_q
!!$  Real(realk), parameter :: PIFAC = 34.986836655249725E0_realk !Two*PI**TwoHalf
!!$  DO iPrimQ = 1, nPrimQ
!!$     q  = expQ(iPrimQ)
!!$     DO iPrimP=1, nPrimP
!!$        p  = expP(iPrimP)
!!$        p_q = p + q
!!$        reducedExponents(iPrimP,iPrimQ) = p*q/p_q
!!$        integralPrefactor(iPrimP,iPrimQ) = PIFAC/(p*q*SQRT(p_q))
!!$     ENDDO
!!$  ENDDO
!!$end subroutine build_reducedExponents_integralPrefactorPQ
    
subroutine build_reducedExponents_integralPrefactorQP(nPrimP,nPrimQ,expQ,expP,&
     & reducedExponents,integralPrefactor)
  implicit none      
  integer,intent(in) :: nPrimP,nPrimQ
  real(realk),intent(in) :: expQ(nPrimQ),expP(nPrimP)
  real(realk),intent(inout) :: reducedExponents(nPrimQ,nPrimP)
  real(realk),intent(inout) :: integralPrefactor(nPrimQ,nPrimP)
  !local variables 
  integer :: iPrimQ,iPrimP
  real(realk) :: p,q,p_q
  Real(realk), parameter :: PIFAC = 34.986836655249725E0_realk !Two*PI**TwoHalf
  DO iPrimP=1, nPrimP
     p  = expP(iPrimP)
     DO iPrimQ = 1, nPrimQ
        q  = expQ(iPrimQ)
        p_q = p + q
        reducedExponents(iPrimQ,iPrimP) = p*q/p_q
        integralPrefactor(iPrimQ,iPrimP) = PIFAC/(p*q*SQRT(p_q))
     ENDDO
  ENDDO
end subroutine build_reducedExponents_integralPrefactorQP

subroutine PrintTypeExpInfo(nPrimP,nPrimQ,reducedExponents,integralPrefactor,lupri)
  implicit none
  integer,intent(in) :: nPrimP,nPrimQ,lupri
  real(realk) :: reducedExponents(nPrimP*nPrimQ),integralPrefactor(nPrimP*nPrimQ)
  !locigal variables
  integer :: iPrimP
  WRITE(lupri,*)'ReducedExponents  QP order'
  do iPrimP=1,nPrimP*nPrimQ
     WRITE(lupri,'(3X,ES18.9)')reducedExponents(iPrimP)
  enddo
  WRITE(lupri,*)'IntegralPrefactor QP order'
  do iPrimP=1,nPrimP*nPrimQ
     WRITE(lupri,'(3X,ES18.9)')integralPrefactor(iPrimP)
  enddo
END subroutine PRINTTYPEEXPINFO
!!$
subroutine PrintTypeInfo(AngmomA,AngmomB,AngmomC,AngmomD,nPrimA,nPrimB,nPrimC,nPrimD,&
     & nContA,nContB,nContC,nContD,expA,ContractCoeffA,expB,ContractCoeffB,&
     & expC,ContractCoeffC,expD,ContractCoeffD,&
     & nAtomsA,nAtomsB,nAtomsC,nAtomsD,ACenter,BCenter,CCenter,DCenter,lupri)
  implicit none
  integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
  integer,intent(in) :: nPrimA,nPrimB,nPrimC,nPrimD,lupri
  integer,intent(in) :: nContA,nContB,nContC,nContD
  integer,intent(in) :: nAtomsA,nAtomsB,nAtomsC,nAtomsD
  real(realk) :: expA(nPrimA),ContractCoeffA(nPrimA,nContA),ACenter(3,nAtomsA)
  real(realk) :: expB(nPrimB),ContractCoeffB(nPrimB,nContB),BCenter(3,nAtomsB)
  real(realk) :: expC(nPrimC),ContractCoeffC(nPrimC,nContC),CCenter(3,nAtomsC)
  real(realk) :: expD(nPrimD),ContractCoeffD(nPrimD,nContD),DCenter(3,nAtomsD)
  !locigal variables
  integer :: iPrimP,iCont,i
  WRITE(lupri,'(2X,A,I1,A,I1,A,I1,A,I1,A)')'Angmom = (',AngmomA,',',AngmomB,',',AngmomC,',',AngmomD,')'
  WRITE(lupri,'(2X,A,I3,A,I3,A,I3,A,I3,A)')'nPrimitives = (',nPrimA,',',nPrimB,',',nPrimC,',',nPrimD,')'
  WRITE(lupri,'(2X,A,I3,A,I3,A,I3,A,I3,A)')'nContracted = (',nContA,',',nContB,',',nContC,',',nContD,')'
  WRITE(lupri,'(2X,A,I1,A,I1,A,I1,A,I1,A)')'Atoms = (',nAtomsA,',',nAtomsB,',',nAtomsC,',',nAtomsD,')'
  WRITE(lupri,'(2X,A)')'ExpA and ContractCoeffA'
  do iPrimP=1,nPrimA
     WRITE(lupri,'(3X,5ES18.9/,(3X,5ES18.9))')expA(iPrimP),(ContractCoeffA(iPrimP,iCont),iCont=1,nContA)
  enddo
  WRITE(lupri,'(2X,A)')'A Centers '
  do i=1,nAtomsA
     WRITE(lupri,'(3ES18.9)') Acenter(1,i),Acenter(2,i),Acenter(3,i)
  enddo

  WRITE(lupri,'(2X,A)')'ExpB and ContractCoeffB'
  do iPrimP=1,nPrimB
     WRITE(lupri,'(3X,5ES18.9/,(3X,5ES18.9))')expB(iPrimP),(ContractCoeffB(iPrimP,iCont),iCont=1,nContB)
  enddo
  WRITE(lupri,'(2X,A)')'B Centers '
  do i=1,nAtomsB
     WRITE(lupri,'(3ES18.9)') Bcenter(1,i),Bcenter(2,i),Bcenter(3,i)
  enddo
  WRITE(lupri,'(2X,A)')'ExpC and ContractCoeffC'
  do iPrimP=1,nPrimC
     WRITE(lupri,'(3X,5ES18.9/,(3X,5ES18.9))')expC(iPrimP),(ContractCoeffC(iPrimP,iCont),iCont=1,nContC)
  enddo
  WRITE(lupri,'(2X,A)')'C Centers '
  do i=1,nAtomsC
     WRITE(lupri,'(3ES18.9)') Ccenter(1,i),Ccenter(2,i),Ccenter(3,i)
  enddo
  WRITE(lupri,'(2X,A)')'ExpD and ContractCoeffD'
  do iPrimP=1,nPrimD
     WRITE(lupri,'(3X,5ES18.9/,(3X,5ES18.9))')expD(iPrimP),(ContractCoeffD(iPrimP,iCont),iCont=1,nContD)
  enddo
  WRITE(lupri,'(2X,A)')'D Centers '
  do i=1,nAtomsD
     WRITE(lupri,'(3ES18.9)') Dcenter(1,i),Dcenter(2,i),Dcenter(3,i)
  enddo
end subroutine PrintTypeInfo

subroutine build_noScreen1(ItypeA,ItypeB,ntypesA,ntypesB,nAtomsA,nAtomsB,nBatchA,nBatchB,&
     & BATCHGAB,MaxGabForTypeCD,THRESHOLD_CS,BatchIndexOfTypeA,BatchIndexOfTypeB,noScreenABin,noScreenABout)
  implicit none
  integer,intent(in) :: ItypeA,ItypeB,nAtomsA,nAtomsB,nBatchA,nBatchB,ntypesA,ntypesB
  real(realk),intent(in) :: BATCHGAB(nBatchA,nBatchB),MaxGabForTypeCD
  real(realk),intent(in) :: THRESHOLD_CS
  integer,intent(in) :: BatchIndexOfTypeA(ntypesA),BatchIndexOfTypeB(ntypesB)
  logical,intent(in) :: noScreenABin(nAtomsA,nAtomsB)
  logical,intent(inout) :: noScreenABout(nAtomsA,nAtomsB)
  !local variables
  integer :: iBatchA,IatomA,iBatchB,iAtomB
  real(realk) :: MAXGAB
!$OMP PARALLEL DO PRIVATE(IatomB,IatomA,iBatchA,iBatchB) FIRSTPRIVATE(nAtomsA,&
!$OMP nAtomsB,MaxGabForTypeCD,THRESHOLD_CS) SHARED(BatchIndexOfTypeA,&
!$OMP BatchIndexOfTypeB,noScreenABout,noScreenABin) SCHEDULE(DYNAMIC,1)
  DO IatomB = 1,nAtomsB
     iBatchB = BatchIndexOfTypeB(ItypeB) + IatomB
     iBatchA = BatchIndexOfTypeA(ItypeA)
     DO IatomA = 1,nAtomsA
        noScreenABout(IatomA,IatomB) = noScreenABin(IatomA,IatomB).AND.&
             &(BATCHGAB(iBatchA+IatomA,iBatchB)*MaxGabForTypeCD.GT.THRESHOLD_CS)
     ENDDO
  ENDDO
!$OMP END PARALLEL DO
end subroutine build_noScreen1

subroutine build_EmptynoScreen1(nAtomsA,nAtomsB,noScreenAB)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB
  logical,intent(inout) :: noScreenAB(nAtomsA,nAtomsB)
  !local variables
  integer :: IatomA,iAtomB
!$OMP PARALLEL DO PRIVATE(IatomB,IatomA) FIRSTPRIVATE(nAtomsA,&
!$OMP nAtomsB) SHARED(noScreenAB) SCHEDULE(DYNAMIC,1)
  DO IatomB = 1,nAtomsB
     DO IatomA = 1,nAtomsA
        noScreenAB(IatomA,IatomB) = .TRUE.
     ENDDO
  ENDDO
!$OMP END PARALLEL DO
end subroutine build_EmptynoScreen1

subroutine copy_noScreen(nAtomsA,nAtomsB,noScreenAB,noScreenAB2)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB
  logical,intent(inout) :: noScreenAB(nAtomsA,nAtomsB)
  logical,intent(inout) :: noScreenAB2(nAtomsA,nAtomsB)
  !local variables
  integer :: IatomA,iAtomB
!$OMP PARALLEL DO PRIVATE(IatomB,IatomA) FIRSTPRIVATE(nAtomsA,&
!$OMP nAtomsB) SHARED(noScreenAB,noScreenAB2) SCHEDULE(DYNAMIC,1)
  DO IatomB = 1,nAtomsB
     DO IatomA = 1,nAtomsA
        noScreenAB2(IatomA,IatomB) = noScreenAB(IatomA,IatomB)
     ENDDO
  ENDDO
!$OMP END PARALLEL DO
end subroutine copy_noScreen


subroutine Build_qcent_Qdistance12_QpreExpFac(nPrimC,nPrimD,nContC,nContD,&
     & expC,expD,Ccenter,Dcenter,ContractCoeffC,ContractCoeffD,Segmented,&
     & qcent,Qdistance12,QpreExpFac,INTPRINT)
  implicit none
  integer,intent(in) :: nPrimC,nPrimD,nContC,nContD,INTPRINT
  real(realk),intent(in) :: expC(nPrimC),expD(nPrimD)
  real(realk),intent(in) :: Ccenter(3),Dcenter(3)
  real(realk),intent(in) :: ContractCoeffC(nPrimC,nContC)
  real(realk),intent(in) :: ContractCoeffD(nPrimD,nContD)
  logical,intent(in) :: Segmented
  real(realk),intent(inout) :: qcent(3,nPrimC,nPrimD),Qdistance12(3)
  real(realk),intent(inout) :: QpreExpFac(nPrimC,nPrimD)
  !local variables
  integer :: i12,i2,i1,offset
  real(realk) :: e2,e1,X,Y,Z,d2,eDX,eDY,eDZ,TMPCCD
  X = Ccenter(1) - Dcenter(1)
  Y = Ccenter(2) - Dcenter(2)
  Z = Ccenter(3) - Dcenter(3)
  d2 = X*X + Y*Y + Z*Z
  Qdistance12(1) = X
  Qdistance12(2) = Y
  Qdistance12(3) = Z
  IF (Segmented) THEN
     DO i2=1,nPrimD
        e2  = expD(i2)       
        eDX = e2*Dcenter(1)
        eDY = e2*Dcenter(2)
        eDZ = e2*Dcenter(3)
        TMPCCD = ContractCoeffD(i2,1)
        DO i1=1,nPrimC
           e1  = expC(i1)
           qcent(1,i1,i2) = (e1*Ccenter(1) + eDX)/(e1+e2)
           Qcent(2,i1,i2) = (e1*Ccenter(2) + eDY)/(e1+e2)
           Qcent(3,i1,i2) = (e1*Ccenter(3) + eDZ)/(e1+e2)
           QpreExpFac(i1,i2) = exp(-e1*e2/(e1+e2)*d2)*ContractCoeffC(i1,1)*TMPCCD
        ENDDO
     ENDDO
  ELSE
     DO i2=1,nPrimD
        e2  = expD(i2)       
        eDX = e2*Dcenter(1)
        eDY = e2*Dcenter(2)
        eDZ = e2*Dcenter(3)
        DO i1=1,nPrimC
           e1  = expC(i1)
           qcent(1,i1,i2) = (e1*Ccenter(1) + eDX)/(e1+e2)
           Qcent(2,i1,i2) = (e1*Ccenter(2) + eDY)/(e1+e2)
           Qcent(3,i1,i2) = (e1*Ccenter(3) + eDZ)/(e1+e2)
           QpreExpFac(i1,i2) = exp(-e1*e2/(e1+e2)*d2)
        ENDDO
     ENDDO
  END IF
end subroutine Build_qcent_Qdistance12_QpreExpFac

subroutine Build_qcent_Qdistance12_QpreExpFacGPU(nPrimC,nPrimD,nContC,nContD,&
     & expC,expD,Ccenter,Dcenter,ContractCoeffC,ContractCoeffD,Segmented,&
     & qcent,Qdistance12,QpreExpFac,INTPRINT)
  implicit none
  integer,intent(in) :: nPrimC,nPrimD,nContC,nContD,INTPRINT
  real(realk),intent(in) :: expC(nPrimC),expD(nPrimD)
  real(realk),intent(in) :: Ccenter(3),Dcenter(3)
  real(realk),intent(in) :: ContractCoeffC(nPrimC,nContC)
  real(realk),intent(in) :: ContractCoeffD(nPrimD,nContD)
  logical,intent(in) :: Segmented
  real(realk),intent(inout) :: qcent(3,nPrimC,nPrimD),Qdistance12(3)
  real(realk),intent(inout) :: QpreExpFac(nPrimC,nPrimD)
  !local variables
  integer :: i12,i2,i1,offset
  real(realk) :: e2,e1,X,Y,Z,d2,eDX,eDY,eDZ,TMPCCD
  IF (Segmented) THEN
!$ACC KERNELS &
!$ACC PRESENT(QpreExpFac,Qcent,Qdistance12,nPrimC,nPrimD,nContC,nContD,&
!$ACC         expC,expD,Ccenter,Dcenter,ContractCoeffC,ContractCoeffD)
     d2 = 0.0E0_realk
     DO I1=1,3
        Qdistance12(I1) = Ccenter(I1) - Dcenter(I1)
        d2 = d2 + Qdistance12(I1)*Qdistance12(I1)
     ENDDO
     DO i2=1,nPrimD
        e2  = expD(i2)       
        eDX = e2*Dcenter(1)
        eDY = e2*Dcenter(2)
        eDZ = e2*Dcenter(3)
        TMPCCD = ContractCoeffD(i2,1)
        DO i1=1,nPrimC
           e1  = expC(i1)
           qcent(1,i1,i2) = (e1*Ccenter(1) + eDX)/(e1+e2)
           Qcent(2,i1,i2) = (e1*Ccenter(2) + eDY)/(e1+e2)
           Qcent(3,i1,i2) = (e1*Ccenter(3) + eDZ)/(e1+e2)
           QpreExpFac(i1,i2) = exp(-e1*e2/(e1+e2)*d2)*ContractCoeffC(i1,1)*TMPCCD
        ENDDO
     ENDDO
!$ACC END KERNELS
  ELSE
!$ACC KERNELS &
!$ACC PRESENT(QpreExpFac,Qcent,Qdistance12,nPrimC,nPrimD,nContC,nContD,&
!$ACC         expC,expD,Ccenter,Dcenter,ContractCoeffC,ContractCoeffD)
     d2 = 0.0E0_realk
     DO I1=1,3
        Qdistance12(I1) = Ccenter(I1) - Dcenter(I1)
        d2 = d2 + Qdistance12(I1)*Qdistance12(I1)
     ENDDO
     DO i2=1,nPrimD
        e2  = expD(i2)       
        eDX = e2*Dcenter(1)
        eDY = e2*Dcenter(2)
        eDZ = e2*Dcenter(3)
        DO i1=1,nPrimC
           e1  = expC(i1)
           qcent(1,i1,i2) = (e1*Ccenter(1) + eDX)/(e1+e2)
           Qcent(2,i1,i2) = (e1*Ccenter(2) + eDY)/(e1+e2)
           Qcent(3,i1,i2) = (e1*Ccenter(3) + eDZ)/(e1+e2)
           QpreExpFac(i1,i2) = exp(-e1*e2/(e1+e2)*d2)
        ENDDO
     ENDDO
!$ACC END KERNELS
  END IF
end subroutine Build_qcent_Qdistance12_QpreExpFacGPU

subroutine Build_Seg_qcent_QpreExpFac(nPrimC,nPrimD,&
     & expC,expD,Ccenter,Dcenter,ContractCoeffC,ContractCoeffD,&
     & qcent,QpreExpFac,INTPRINT)
  implicit none
  integer,intent(in) :: nPrimC,nPrimD,INTPRINT
  real(realk),intent(in) :: expC(nPrimC),expD(nPrimD)
  real(realk),intent(in) :: Ccenter(3),Dcenter(3)
  real(realk),intent(in) :: ContractCoeffC(nPrimC)
  real(realk),intent(in) :: ContractCoeffD(nPrimD)
  real(realk),intent(inout) :: qcent(3,nPrimC,nPrimD)
  real(realk),intent(inout) :: QpreExpFac(nPrimC,nPrimD)
  !local variables
  integer :: i12,i2,i1,offset
  real(realk) :: e2,e1,X,Y,Z,d2,eDX,eDY,eDZ,TMPCCD
  X = Ccenter(1) - Dcenter(1)
  Y = Ccenter(2) - Dcenter(2)
  Z = Ccenter(3) - Dcenter(3)
  d2 = X*X + Y*Y + Z*Z
  DO i2=1,nPrimD
     e2  = expD(i2)       
     eDX = e2*Dcenter(1)
     eDY = e2*Dcenter(2)
     eDZ = e2*Dcenter(3)
     TMPCCD = ContractCoeffD(i2)
     DO i1=1,nPrimC
        e1  = expC(i1)
        qcent(1,i1,i2) = (e1*Ccenter(1) + eDX)/(e1+e2)
        Qcent(2,i1,i2) = (e1*Ccenter(2) + eDY)/(e1+e2)
        Qcent(3,i1,i2) = (e1*Ccenter(3) + eDZ)/(e1+e2)
        QpreExpFac(i1,i2) = exp(-e1*e2/(e1+e2)*d2)*ContractCoeffC(i1)*TMPCCD
     ENDDO
  ENDDO
end subroutine Build_Seg_qcent_QpreExpFac

SUBROUTINE Build_pcent_Pdistance12_PpreExpFac(nPrimA,nPrimB,natomsA,natomsB,nContA,nContB,&
     & inversexpP,expA,expB,Acenter,Bcenter,ContractCoeffA,ContractCoeffB,Segmented,&
     & pcentPass,Pdistance12Pass,PpreExpFacPass,INTPRINT)
  implicit none
  integer,intent(in) :: nPrimA,nPrimB,natomsA,natomsB,nContA,nContB,INTPRINT
  real(realk),intent(in) :: inversexpP(nPrimA,nPrimB),expA(nPrimA),expB(nPrimB)
  real(realk),intent(in) :: Acenter(3,natomsA),Bcenter(3,natomsB)
  real(realk),intent(in) :: ContractCoeffA(nPrimA,nContA)
  real(realk),intent(in) :: ContractCoeffB(nPrimB,nContB)
  logical,intent(in)     :: Segmented!,noScreenAB(nAtomsA,nAtomsB)
  real(realk),intent(inout) :: PcentPass(3,nPrimA,nPrimB,natomsA,natomsB)
  real(realk),intent(inout) :: Pdistance12Pass(3,natomsA,natomsB)
  real(realk),intent(inout) :: PpreExpFacPass(nPrimA,nPrimB,natomsA,natomsB)
  !local variables
  integer :: i12,i2,i1,offset,IatomA,IatomB
  real(realk) :: e2,e1,X,Y,Z,d2,AX,AY,AZ,BX,BY,BZ,TMPCCB,tmpe2d2,eBX,eBY,eBZ
  IF (Segmented) THEN
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(IatomB,BX,BY,BZ,IatomA,AX,AY,AZ,X,Y,Z,d2,&
!$OMP e2,e1,eBX,eBY,eBZ,tmpe2d2,TMPCCB,i1,i2) FIRSTPRIVATE(nAtomsB,nAtomsA,nPrimA,&
!$OMP nPrimB) SHARED(expA,expB,inversexpP,Pdistance12Pass,&
!$OMP Acenter,Bcenter,pcentPass,ContractCoeffA,ContractCoeffB,PpreExpFacPass) SCHEDULE(DYNAMIC,1)
  DO IatomB = 1,nAtomsB
   BX = Bcenter(1,IatomB)
   BY = Bcenter(2,IatomB)
   BZ = Bcenter(3,IatomB)
   DO IatomA = 1,nAtomsA
     AX = Acenter(1,IatomA)
     AY = Acenter(2,IatomA)
     AZ = Acenter(3,IatomA)
     X = AX - BX
     Y = AY - BY
     Z = AZ - BZ
     Pdistance12Pass(1,iAtomA,IatomB) = X
     Pdistance12Pass(2,iAtomA,IatomB) = Y
     Pdistance12Pass(3,iAtomA,IatomB) = Z
     d2 = X*X + Y*Y + Z*Z
     DO i2=1,nPrimB
      e2  = expB(i2)       
      eBX = e2*BX
      eBY = e2*BY
      eBZ = e2*BZ
      tmpe2d2 = e2*d2
      TMPCCB = ContractCoeffB(i2,1)
      DO i1=1,nPrimA
        pcentPass(1,i1,i2,iAtomA,IatomB) = (AX*expA(i1) + eBX)*inversexpP(i1,i2)
        pcentPass(2,i1,i2,iAtomA,IatomB) = (AY*expA(i1) + eBY)*inversexpP(i1,i2)
        pcentPass(3,i1,i2,iAtomA,IatomB) = (AZ*expA(i1) + eBZ)*inversexpP(i1,i2)
        PpreExpFacPass(i1,i2,iAtomA,IatomB) = exp(-expA(i1)*tmpe2d2*inversexpP(i1,i2))*ContractCoeffA(i1,1)*TMPCCB
      ENDDO
     ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO 
ELSE
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(IatomB,BX,BY,BZ,IatomA,AX,AY,AZ,X,Y,Z,d2,&
!$OMP e2,e1,eBX,eBY,eBZ,tmpe2d2,TMPCCB,i1,i2) FIRSTPRIVATE(nAtomsB,nAtomsA,nPrimA,&
!$OMP nPrimB) SHARED(expA,expB,inversexpP,Pdistance12Pass,&
!$OMP Acenter,Bcenter,pcentPass,ContractCoeffA,PpreExpFacPass) SCHEDULE(DYNAMIC,1)
  DO IatomB = 1,nAtomsB
   BX = Bcenter(1,IatomB)
   BY = Bcenter(2,IatomB)
   BZ = Bcenter(3,IatomB)
   DO IatomA = 1,nAtomsA
     AX = Acenter(1,IatomA)
     AY = Acenter(2,IatomA)
     AZ = Acenter(3,IatomA)
     X = AX - BX
     Y = AY - BY
     Z = AZ - BZ
     Pdistance12Pass(1,iAtomA,IatomB) = X
     Pdistance12Pass(2,iAtomA,IatomB) = Y
     Pdistance12Pass(3,iAtomA,IatomB) = Z
     d2 = X*X + Y*Y + Z*Z
     DO i2=1,nPrimB
      e2  = expB(i2)       
      eBX = e2*BX
      eBY = e2*BY
      eBZ = e2*BZ
      tmpe2d2 = e2*d2
      DO i1=1,nPrimA
       pcentPass(1,i1,i2,iAtomA,IatomB) = (AX*expA(i1) + eBX)*inversexpP(i1,i2)
       pcentPass(2,i1,i2,iAtomA,IatomB) = (AY*expA(i1) + eBY)*inversexpP(i1,i2)
       pcentPass(3,i1,i2,iAtomA,IatomB) = (AZ*expA(i1) + eBZ)*inversexpP(i1,i2)
       PpreExpFacPass(i1,i2,iAtomA,IatomB) = exp(-tmpe2d2*expA(i1)*inversexpP(i1,i2))
      ENDDO
     ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO 
ENDIF
end SUBROUTINE Build_pcent_Pdistance12_PpreExpFac

SUBROUTINE Build_pcent_PpreExpFac(nPrimA,nPrimB,natomsA,natomsB,nContA,nContB,&
     & inversexpP,expA,expB,Acenter,Bcenter,ContractCoeffA,ContractCoeffB,Segmented,&
     & pcentPass,PpreExpFacPass,INTPRINT)
  implicit none
  integer,intent(in) :: nPrimA,nPrimB,natomsA,natomsB,nContA,nContB,INTPRINT
  real(realk),intent(in) :: inversexpP(nPrimA,nPrimB),expA(nPrimA),expB(nPrimB)
  real(realk),intent(in) :: Acenter(3,natomsA),Bcenter(3,natomsB)
  real(realk),intent(in) :: ContractCoeffA(nPrimA,nContA)
  real(realk),intent(in) :: ContractCoeffB(nPrimB,nContB)
  logical,intent(in)     :: Segmented!,noScreenAB(nAtomsA,nAtomsB)
  real(realk),intent(inout) :: PcentPass(3,nPrimA,nPrimB,natomsA,natomsB)
  real(realk),intent(inout) :: PpreExpFacPass(nPrimA,nPrimB,natomsA,natomsB)
  !local variables
  integer :: i12,i2,i1,offset,IatomA,IatomB
  real(realk) :: e2,e1,X,Y,Z,d2,AX,AY,AZ,BX,BY,BZ,TMPCCB,tmpe2d2,eBX,eBY,eBZ
  IF (Segmented) THEN
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(IatomB,BX,BY,BZ,IatomA,AX,AY,AZ,X,Y,Z,d2,&
!$OMP e2,e1,eBX,eBY,eBZ,tmpe2d2,TMPCCB,i1,i2) FIRSTPRIVATE(nAtomsB,nAtomsA,nPrimA,&
!$OMP nPrimB) SHARED(expA,expB,inversexpP,&
!$OMP Acenter,Bcenter,pcentPass,ContractCoeffA,ContractCoeffB,PpreExpFacPass) SCHEDULE(DYNAMIC,1)
  DO IatomB = 1,nAtomsB
   BX = Bcenter(1,IatomB)
   BY = Bcenter(2,IatomB)
   BZ = Bcenter(3,IatomB)
   DO IatomA = 1,nAtomsA
     AX = Acenter(1,IatomA)
     AY = Acenter(2,IatomA)
     AZ = Acenter(3,IatomA)
     X = AX - BX
     Y = AY - BY
     Z = AZ - BZ
     d2 = X*X + Y*Y + Z*Z
     DO i2=1,nPrimB
      e2  = expB(i2)       
      eBX = e2*BX
      eBY = e2*BY
      eBZ = e2*BZ
      tmpe2d2 = e2*d2
      TMPCCB = ContractCoeffB(i2,1)
      DO i1=1,nPrimA
        pcentPass(1,i1,i2,iAtomA,IatomB) = (AX*expA(i1) + eBX)*inversexpP(i1,i2)
        pcentPass(2,i1,i2,iAtomA,IatomB) = (AY*expA(i1) + eBY)*inversexpP(i1,i2)
        pcentPass(3,i1,i2,iAtomA,IatomB) = (AZ*expA(i1) + eBZ)*inversexpP(i1,i2)
        PpreExpFacPass(i1,i2,iAtomA,IatomB) = exp(-expA(i1)*tmpe2d2*inversexpP(i1,i2))*ContractCoeffA(i1,1)*TMPCCB
      ENDDO
     ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO 
ELSE
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(IatomB,BX,BY,BZ,IatomA,AX,AY,AZ,X,Y,Z,d2,&
!$OMP e2,e1,eBX,eBY,eBZ,tmpe2d2,TMPCCB,i1,i2) FIRSTPRIVATE(nAtomsB,nAtomsA,nPrimA,&
!$OMP nPrimB) SHARED(expA,expB,inversexpP,&
!$OMP Acenter,Bcenter,pcentPass,ContractCoeffA,PpreExpFacPass) SCHEDULE(DYNAMIC,1)
  DO IatomB = 1,nAtomsB
   BX = Bcenter(1,IatomB)
   BY = Bcenter(2,IatomB)
   BZ = Bcenter(3,IatomB)
   DO IatomA = 1,nAtomsA
     AX = Acenter(1,IatomA)
     AY = Acenter(2,IatomA)
     AZ = Acenter(3,IatomA)
     X = AX - BX
     Y = AY - BY
     Z = AZ - BZ
     d2 = X*X + Y*Y + Z*Z
     DO i2=1,nPrimB
      e2  = expB(i2)       
      eBX = e2*BX
      eBY = e2*BY
      eBZ = e2*BZ
      tmpe2d2 = e2*d2
      DO i1=1,nPrimA
       pcentPass(1,i1,i2,iAtomA,IatomB) = (AX*expA(i1) + eBX)*inversexpP(i1,i2)
       pcentPass(2,i1,i2,iAtomA,IatomB) = (AY*expA(i1) + eBY)*inversexpP(i1,i2)
       pcentPass(3,i1,i2,iAtomA,IatomB) = (AZ*expA(i1) + eBZ)*inversexpP(i1,i2)
       PpreExpFacPass(i1,i2,iAtomA,IatomB) = exp(-tmpe2d2*expA(i1)*inversexpP(i1,i2))
      ENDDO
     ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO 
ENDIF
end SUBROUTINE Build_pcent_PpreExpFac

SUBROUTINE Build_pcent_Pdistance12_PpreExpFac2(nPrimP,nPasses,&
     & PcentPass,Pdistance12Pass,PpreExpFacPass,&
     & pcent,Pdistance12,PpreExpFac,nAtomsA,nAtomsB,IatomA,IatomB)
  implicit none
  integer,intent(in) :: nPrimP,nPasses,nAtomsA,nAtomsB,IatomA,IatomB
  real(realk),intent(in) :: PcentPass(3,nPrimP,nAtomsA,nAtomsB),Pdistance12Pass(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: PpreExpFacPass(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: Pcent(3,nPrimP),Pdistance12(3)
  real(realk),intent(inout) :: PpreExpFac(nPrimP)
  !local variables
  integer :: i
  Pdistance12(1) = Pdistance12Pass(1,iAtomA,iAtomB)
  Pdistance12(2) = Pdistance12Pass(2,iAtomA,iAtomB)
  Pdistance12(3) = Pdistance12Pass(3,iAtomA,iAtomB)
  DO i=1,nPrimP
     pcent(1,i) = pcentPass(1,i,iAtomA,iAtomB)
     pcent(2,i) = pcentPass(2,i,iAtomA,iAtomB)
     pcent(3,i) = pcentPass(3,i,iAtomA,iAtomB)
     PpreExpFac(i) = PpreExpFacPass(i,iAtomA,iAtomB) 
  ENDDO
end SUBROUTINE Build_pcent_Pdistance12_PpreExpFac2

subroutine build_ichor_AOextent(MaxnAtomsA,MaxnprimA,MaxnContA,ntypesA,exponentsOfTypeA,ODscreen,&
     & nContOfTypeA,nPrimOfTypeA,ContractCoeffOfTypeA,AngmomOfTypeA,OD_Threshold,ExtentOfTypeA)
  implicit none
  integer,intent(in) :: nTypesA,MaxnAtomsA,MaxnPrimA,MaxnContA
  Integer,intent(in) :: AngmomOfTypeA(ntypesA),nContOfTypeA(ntypesA),nPrimOfTypeA(ntypesA)
  Real(realk),intent(in) :: exponentsOfTypeA(MaxnprimA,ntypesA)
  Real(realk),intent(in) :: ContractCoeffOfTypeA(MaxnprimA,MaxnContA,ntypesA)
  Real(realk),intent(in) :: OD_Threshold
  Real(realk),intent(inout) :: ExtentOfTypeA(ntypesA)
  !local variables
  real(realk) :: extent,extent2,maxContraction,r2,fun,funD,rold,rnew,expo
  integer :: nPrim,Angmom,nCont,iPrim,iCont,itypeA
  logical :: ODscreen
  IF(ODscreen)THEN
   DO ItypeA=1,nTypesA     
    extent2 = 0E0_realk
    Angmom = AngmomOfTypeA(ItypeA)
    nPrim = nPrimOfTypeA(ItypeA)
    nCont = nContOfTypeA(ItypeA)
    DO iPrim=1,nPrim
     maxContraction = tiny(1E0_realk)
     DO iCont=1,nCont
        maxContraction = max(maxContraction,abs(ContractCoeffOfTypeA(iPrim,iCont,itypeA)))
     ENDDO
     expo = exponentsOfTypeA(iPrim,itypeA)
     IF((OD_Threshold.LT.tiny(1E0_realk)).OR.(nCont.LT. 1).OR.(expo.LT.tiny(1E0_realk))) THEN
        WRITE(*,*) 'Error in build_ichor_AOextent', OD_Threshold, nCont,expo
        CALL ICHORQUIT('Error in build_ichor_AOextent',-1)
     ENDIF
     !r2 = (-log(Threshold)+log(maxContraction)+log(real(nCont)))/Exponents(iPrim)
     r2 = (-log(OD_Threshold)+log(maxContraction))/Expo
     !   Take the above expression to be a constant A. We should then solve 
     !      r^2 = A + l/a*ln(r)       (i)
     !   to get a proper extent. We instead take A to be an estimate of r^2
     !   and make the correction r^2 = A + l/a*ln(sqrt(A)). The equation (i)
     !   can of course instead be solved iteratively.
     IF (r2.GT. 0E0_realk) THEN
      IF(angmom .GT. 0)THEN
       Rold = sqrt(r2 + angmom*log(sqrt(r2))/Expo)
       fun = ABS(maxContraction*(rold**angmom)*exp(-Expo*Rold**2))
       IF(OD_threshold .LE. fun)THEN
        DO
         funD = ABS( maxContraction*angmom*(rold**(angmom-1))*&
              &exp(-Expo*Rold**2)+maxContraction*(rold**angmom)&
              &*(-2*Expo*Rold)*exp(-Expo*Rold**2))
         Rnew = ROLD - (OD_Threshold-fun)/funD
         IF(ABS(Rnew-ROLD) .LE. 1.0E-24_8)THEN
            write(*,*)'calculation stalled in build_ichor_AOextent'
            EXIT
         ENDIF
         fun = ABS(maxContraction*(Rnew)**angmom*exp(-Expo*Rnew**2))
         Rold=Rnew
         IF(ABS(fun-OD_Threshold).LE. OD_Threshold*1.0E-11_8)EXIT
         IF(Rold.NE.Rold)THEN
            write(*,*)'build_ichor_AOextent: found NaN in aobatch iteratvie extent'
            EXIT
         ENDIF
        ENDDO
       ENDIF
       r2=Rold*Rold
      ENDIF
     ENDIF
     extent2 = max(extent2,r2)
    ENDDO
    IF (extent2.LT. 0E0_realk) THEN
     WRITE(*,*) 'Negative squared distance in build_ichor_AOextent',extent2
     CALL ICHORQUIT('Negative squared distance in build_ichor_AOextent',-1)
    ENDIF
    ExtentOfTypeA(itypeA) = sqrt(extent2)
   enddo
  ELSE
   DO ItypeA=1,nTypesA     
    ExtentOfTypeA(itypeA) = 0E0_realk
   ENDDO
  ENDIF
end subroutine build_ichor_AOextent

subroutine ODscreen_noScreen(nAtomsC,nAtomsD,Ccenter,Dcenter,&
     & extent2CD,noScreenCD)
  implicit none
  integer,intent(in) :: nAtomsC,nAtomsD
  real(realk),intent(in) :: Ccenter(3,nAtomsC),Dcenter(3,nAtomsD)
  real(realk),intent(in) :: extent2CD
  logical,intent(inout)  :: noScreenCD(nAtomsC,nAtomsD)
  !local variables
  integer :: IatomD,IatomC
  real(realk) :: DX,DY,DZ,X,Y,Z,distance2
!$OMP PARALLEL DO PRIVATE(IatomD,IatomC,DX,DY,DZ,X,Y,Z,&
!$OMP distance2) SHARED(noScreenCD) FIRSTPRIVATE(nAtomsD,&
!$OMP nAtomsC,extent2CD) SCHEDULE(DYNAMIC,1)
  DO IatomD = 1,nAtomsD
     DX = -Dcenter(1,IatomD)
     DY = -Dcenter(2,IatomD)
     DZ = -Dcenter(3,IatomD)
     DO IatomC = 1,nAtomsC
        X = Ccenter(1,IatomC) + DX
        Y = Ccenter(2,IatomC) + DY
        Z = Ccenter(3,IatomC) + DZ
        distance2 = X*X + Y*Y + Z*Z
        noScreenCD(iAtomC,iAtomD) = distance2.LE.extent2CD
     ENDDO
  ENDDO
!$OMP END PARALLEL DO 
END subroutine ODscreen_noScreen

subroutine ichorzero(dx, length)
  implicit none
  !Length of array
  integer, intent(in)      :: length
  !Array to be nullified
  real(realk), intent(inout) :: dx(length)
  !local
  integer                  :: i
!!$OMP PARALLEL DO PRIVATE(I) FIRSTPRIVATE(length) SHARED(dx) SCHEDULE(DYNAMIC,13)
  do i = 1, length
     dx(i) = 0.0E0_realk
  enddo
!!$OMP END PARALLEL DO
end subroutine ichorzero

subroutine ichorzero5(OutputStorage, Dim1,Dim2,Dim3,Dim4,Dim5)
  implicit none
  !Length of array
  Integer,intent(in) :: Dim1,Dim2,Dim3,Dim4,Dim5
  !Array to be nullified
  real(realk), intent(inout) :: OutputStorage(Dim1,Dim2,Dim3,Dim4,Dim5)
  !local
  integer                  :: i,j,k,l,m
!$OMP PARALLEL DEFAULT(none) PRIVATE(I,J,k,l,m) FIRSTPRIVATE(dim1,&
!$OMP dim2,dim3,dim4,dim5) SHARED(OutputStorage)
  do m = 1, dim5
   do l = 1, dim4
    do k = 1, dim3
!$OMP DO SCHEDULE(DYNAMIC,3)
     do j = 1, dim2
      do i = 1, dim1
       OutputStorage(i,j,k,l,m) = 0.0E0_realk
      enddo
     enddo
!$OMP END DO NOWAIT
    enddo
   enddo
  enddo
!$OMP END PARALLEL
end subroutine ichorzero5

subroutine ichorzero2(OutputStorage, Dim1,Dim2)
  implicit none
  !Length of array
  Integer,intent(in) :: Dim1,Dim2
  !Array to be nullified
  real(realk), intent(inout) :: OutputStorage(Dim1,Dim2)
  !local
  integer :: i,j
  logical :: moda,modb
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(I,J) SHARED(OutputStorage) FIRSTPRIVATE(dim1,&
!$OMP dim2) SCHEDULE(DYNAMIC,13)
  do j=1,dim2
     do i=1,dim1
        OutputStorage(i,j)=0.0E0_realk
     enddo
  enddo
!$OMP END PARALLEL DO 
end subroutine ichorzero2

END MODULE IchorCommonModule
