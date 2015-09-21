module files

    logical, save :: first = .true.
    logical, save :: access_stream = .false.
    integer, save :: IUNTAB(99)

public :: lsopen,lsclose,access_stream
private
contains 

!> \brief Routine for opening files in the DALTON program.
!> \author S. Host
!> \date June 2010
!>
!> Based on LSOPEN by K. Ruud
!> General purpose routine for opening files in the Dalton program.
!> The routine will dynamically allocate unit numbers that will
!> become available again when the file is closed using lsclose.
!>
!> Direct access files are no longer supported!
!>
!> These files are strongly machine dependent, although care has been 
!> taken to avoid using unit numbers that are illegal or reserved on a
!> particular architecture. However, the use of lsopen and lsclose should
!> remove much of the machine dependence in the rest of the Dalton 
!> program
!>
!> Input:
!>  LUNIT    Suggested unit number (OPTIONAL, but is mandatory if an unnamed
!>           file is reopened after having been closed with STATUS='KEEP')
!>           Otherwise it is recommended to not assign this.
!>  FILEIN   Suggested name for the file (OPTIONAL, but strongly recommended)
!>  STATIN   Suggested status of the file (OPTIONAL and maybe not recommended)
!>  FORMIN   Formatted or unformatted file format. Default is 'UNFORMATTED'
!>  POSIN    File positioning (OPTIONAL, default is 'ASIS')
!>
!> Output:
!>  LUNIT    Assigned file unit number
!>
subroutine lsopen(lunit,filein,statin,formin,posin)
   implicit none
   !> Suggested unit number
   integer, intent(inout)  :: lunit
   !> Filename
   character(len=*), intent(in) :: filein
   !> Status of the file ('NEW','OLD','SCRATCH','UNKNOWN'(default))
   character(len=*), intent(in) :: statin 
   !> Formatted or unformatted file format ('FORMATTED','UNFORMATTED'(default))
   character(len=*), intent(in) :: formin 
   !> Optional file positioning ('ASIS'(default),'REWIND','APPEND')
   character(len=*), optional, intent(in) :: posin
   integer       :: FILELEN, STATLEN, FORMLEN, iunit
   integer       :: LENOUT,LENWRK,IOS,i
   character(len=200) :: filename, filestatus, fileformat 
   character(len=20) :: outfil,acc_type,pos
   character(len=200) :: wrkdir
   logical :: fileexists

   ! default position is 'ASIS'
   pos='ASIS'
   if (present(posin)) pos = posin

   ! set access type
   if (access_stream) then
      acc_type='STREAM'
   else
      acc_type='SEQUENTIAL'
   end if

   do I=1,200
      filestatus(I:I) = ' '
      filename(I:I) = ' '
      fileformat(I:I) = ' '
   enddo
   if (first) then
      call ls_izero(IUNTAB,99)
#if defined(SYS_HPUX)
      !
      !     This avoids us writing to unit 7 of HP-systems (which corresponds to
      !     standard error on this machine)
      !
      IUNTAB(7) = 1
#endif
      first = .false.
   endif
   !
   !     We first deal with the unit number
   !
   FILELEN = LEN(FILEIN)
   STATLEN = LEN(STATIN)
   FORMLEN = LEN(FORMIN)
   filename(1:FILELEN) = FILEIN(1:FILELEN)
   filestatus(1:STATLEN) = STATIN(1:STATLEN)
   fileformat(1:FORMLEN) = FORMIN(1:FORMLEN)

   if ((LUNIT < 1) .or. (LUNIT > 99)) then
      !
      !     Unit number left unassigned, we get to decide!
      !
      !Quickfix: Start at unit number 50, since some of the first numbers
      !are used by "old" dalton. Not really pretty, but can be removed
      !once we are completely free of the old code.
      iunit = 50
      do
         iunit = iunit + 1
         IF (iunit > 99) call error_open(9001,lunit,filename,filelen)
         IF ((iunit == 5) .or. (iunit == 6)) then
            cycle
         endif
         IF (IUNTAB(iunit) /= 0) then
            cycle
         endif
         LUNIT = iunit
         exit
      enddo
   else
      !
      !     The user has requested a specific unit number. We don't 
      !     quite trust the user, so we check that
      !     1) It is not unit 5, or 6
      !     2) The file either has been closed with status='KEEP' or is
      !        not currently in use
      !
      IF ((LUNIT == 5) .or. (LUNIT == 6)) call error_open(9002,lunit,filename,filelen)
      IF (IUNTAB(LUNIT) == 1) call error_open(9003,lunit,filename,filelen)
   END IF

   IF (filename(1:FILELEN) == 'LSDALTON.OUT') THEN
      OUTFIL = '                    '
#if defined(SYS_T3D)||defined(SYS_T90)
      CALL PXFGETENV ('OUTFIL',6,OUTFIL,LOUTFL,IERR)
#else 
      CALL GETENV ('OUTFIL',OUTFIL)
#endif
      IF (OUTFIL(1:1) /= ' ') THEN
         LENOUT = 0
         DO I = 1, 20
            IF (OUTFIL(I:I) .EQ. ' ') exit
            LENOUT = LENOUT + 1
         END DO

#if defined(SYS_T3D)||defined(SYS_T90)
         CALL PXFGETENV ('WRKDIR',6,WRKDIR,LWRKDR,IERR)
#else 
         CALL GETENV ('WRKDIR',WRKDIR)
#endif
         LENWRK = 0
         DO I = 1, 200
            IF (WRKDIR(I:I) .EQ. ' ') exit
            LENWRK = LENWRK + 1
         END DO

         FILELEN = LENWRK + LENOUT + 1
         filename(1:FILELEN) = WRKDIR(1:LENWRK)//'/'//OUTFIL(1:LENOUT)
      END IF
   END IF
   IF (STATLEN .EQ. 1) THEN
      STATLEN = 7
      filestatus(1:7) = 'UNKNOWN'
   END IF

   !
   !     We've got a file number now
   !
   IF (filestatus(1:3) == 'OLD' .and. IUNTAB(LUNIT) == 0) THEN
      !
      !        This better be a file with a name, and it better exist
      !
      IF (FILELEN == 1) call error_open(9005,lunit,filename,filelen)

      INQUIRE(FILE=filename(1:FILELEN),EXIST=fileexists,IOSTAT=IOS)

      IF (.NOT. fileexists) call error_open(9007,lunit,filename,filelen)

      IF (FORMLEN == 1) THEN
         FORMLEN = 11
         fileformat(1:FORMLEN)='UNFORMATTED'
      END IF

      OPEN(UNIT=LUNIT,FILE=filename(1:FILELEN),STATUS='OLD', &
         &     FORM=fileformat(1:FORMLEN),IOSTAT=IOS,ACCESS=ACC_TYPE,POSITION=POS)

      if (ios /= 0) call error_open(9004,lunit,filename,filelen,ios)

   ELSE
      IF (filestatus(1:3) == 'NEW' .AND. IUNTAB(LUNIT) /= 0) call error_open(9006,lunit,filename,filelen)

      if (filename == ' ') then
         call lsquit('LSOPEN: file must have a name!',-1)
      endif
      !         IF (FILELEN == 1) THEN
      !            FILELEN = 6
      !            FILENM = 'UNIT'//CHRNOS(LUNIT/10)//
      !     &               CHRNOS(MOD(LUNIT,10))
      !         END IF

      IF (FORMLEN == 1) THEN
         FORMLEN = 11
         fileformat(1:11)='UNFORMATTED'
      END IF

      IF (filestatus(1:7) == 'SCRATCH') THEN 
         OPEN(UNIT=LUNIT,STATUS='SCRATCH',FORM=fileformat(1:FORMLEN),IOSTAT=IOS, &
            & ACCESS=ACC_TYPE, POSITION=POS)
         if (ios /= 0) call error_open(9004,lunit,filename,filelen,ios)
      ELSE

         IF(filestatus(1:STATLEN) == 'NEW') THEN
            !hjaaj/may2000:... if filename already exists the program will abort
            !                  thus we inquire first and delete if necessary.
            !                  This will often be the case if we restart a
            !                  calculation.

            INQUIRE(FILE=filename(1:FILELEN),EXIST=fileexists,IOSTAT=IOS)

            IF (fileexists) THEN
               OPEN(UNIT=LUNIT,FILE=filename(1:FILELEN),STATUS='OLD', &
                  & FORM=fileformat(1:FORMLEN),IOSTAT=IOS,ACCESS=ACC_TYPE,POSITION=POS)
               if (ios /= 0) call error_open(9004,lunit,filename,filelen,ios)
               CLOSE(UNIT=LUNIT,STATUS='DELETE')
            END IF
         END IF

         OPEN(UNIT=LUNIT,FILE=filename(1:FILELEN),STATUS=filestatus(1:STATLEN), &
            & FORM=fileformat(1:FORMLEN),IOSTAT=IOS,ACCESS=ACC_TYPE,POSITION=POS)
         if (ios /= 0) call error_open(9004,lunit,filename,filelen,ios)
      END IF
   ENDIF

   IUNTAB(LUNIT) = 1

end subroutine lsopen

subroutine error_open(error,lunit,filename,filelen,ios)
implicit none
   integer, intent(in) :: error, lunit, filelen
   character(len=80),intent(in) :: filename
   integer,intent(in),optional  :: ios

   if (error == 9001) then
      WRITE (6,'(//A/A/A//A)') &
     &   '--> ERROR (LSOPEN) NO MORE AVAILABLE FILENUMBERS!', &
     &   '--> THIS CALCULATION EITHER NEEDS TOO MANY SIMULTANEOUS '// &
     &   'FILES OR', &
     &   '--> SOMEBODY HAS FORGOTTEN TO CLOSE FILES IN THE SOURCE CODE', &
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
      !CALL QTRACE(6)
      CALL lsQUIT('ERROR (LSOPEN) NO MORE FILE NUMBERS',-1)
   else if (error == 9002) then
      WRITE (6,'(//A/A,I3/A//2A//A)') &
     &   '--> ERROR (LSOPEN) TRYING TO OPEN AN ILLEGAL FILE NUMBER', &
     &   '--> SOMEBODY HAS TRIED TO OPEN UNIT',LUNIT, &
     &   '--> THE PROGRAM DOES NOT ALLOW THE USE OF THIS RESERVED'// &
     &   'UNIT NUMBER', &
     &   '--> Name of offending file (if any): ',filename(1:FILELEN), &
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
      !CALL QTRACE(6)
      CALL lsQUIT('ERROR (LSOPEN) ILLEGAL FILE NUMBER REQUESTED',-1)
   else if (error == 9003) then
      WRITE (6,'(//A/A//A,I5/2A//A)') &
     &   '--> ERROR (LSOPEN) TRYING TO OPEN A FILE ALREADY IN USE', &
     &   '--> SOMEBODY IS TRYING TO USE A FILENUMBER THAT IS '// &
     &   'ALREADY IN USE', &
     &   '--> Offending UNIT number: ',LUNIT, &
     &   '--> Name of offending file (if any): ', filename(1:FILELEN), &
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
      !CALL QTRACE(6)
      CALL lsQUIT('ERROR (LSOPEN) TRYING TO OPEN A FILE ALREADY IN USE',-1)
   else if (error == 9004) then
     if (present(ios)) then
         WRITE (6,'(//A,I3/2A/A,I7)') &
        &   '--> ERROR (LSOPEN) UPON TRYING TO OPEN FILE ON UNIT',LUNIT, &
        &   '--> with filename ',filename(1:FILELEN), &
        &   '--> IOSTAT ERROR CODE RETURNED ',IOS
         !CALL QTRACE(6)
         CALL lsQUIT('ERROR (LSOPEN) UPON OPENING A FILE',-1)
     else
        CALL lsQUIT('ios must be present when error=9004',-1)
     endif
   else if (error == 9005) then
      WRITE (6,'(//A/A/A/A,I5/A,I5/2A//A)') &
     &   '--> ERROR (LSOPEN) TRYING TO OPEN A NON-EXISTING OLD FILE', &
     &   '--> A FILE HAS BEEN SPECIFIED TO BE OLD, BUT THE TABLE', &
     &   '--> ENTRY INDICATES THAT IT DOES NOT EXIST', &
     &   '--> IOSTAT ERROR CODE RETURNED ',IOS,&
     &   '--> Offending UNIT number: ',LUNIT, &
     &   '--> Name of offending file (if any): ',filename(1:FILELEN), &
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
      !CALL QTRACE(6)
      CALL lsQUIT('ERROR (LSOPEN) TRYING TO OPEN A NON-EXISTING FILE AS OLD',-1)
   else if (error == 9006) then
      WRITE (6,'(//A/A/A//A,I5/2A//A)') &
     &   '--> ERROR (LSOPEN) TRYING TO OPEN AN EXISTING NEW FILE', &
     &   '--> A FILE HAS BEEN SPECIFIED TO BE NEW, BUT THE TABLE', &
     &   '--> ENTRY INDICATES THAT IT ALREADY EXISTS', &
     &   '--> Offending UNIT number: ',LUNIT, &
     &   '--> Name of offending file (if any): ',filename(1:FILELEN), &
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
      !CALL QTRACE(6)
      CALL lsQUIT('ERROR (LSOPEN) TRYING TO OPEN AN EXISTING FILE AS NEW',-1)
   else if (error == 9007) then
      print*,'IOSTAT ERROR CODE RETURNED ',IOS
      print*,'--> Offending UNIT number: ',LUNIT
      print*,'--> Name of offending file (if any): ',filename(1:FILELEN)
      WRITE (6,'(//A/A/A/A,I5/A,I5/2A//A)') &
     &   '--> ERROR (LSOPEN) TRYING TO OPEN A NON-EXISTING OLD FILE', &
     &   '--> A FILE HAS BEEN SPECIFIED TO BE OLD,', &
     &   '--> BUT THE FILE DOES NOT EXIST', &
     &   '--> IOSTAT ERROR CODE RETURNED ',IOS,&
     &   '--> Offending UNIT number: ',LUNIT, &
     &   '--> Name of offending file (if any): ',filename(1:FILELEN), &
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
     print*,'--> ERROR (LSOPEN) TRYING TO OPEN A NON-EXISTING OLD FILE', &
     &   '--> A FILE HAS BEEN SPECIFIED TO BE OLD,', &
     &   '--> BUT THE FILE DOES NOT EXIST', &
     &   '--> IOSTAT ERROR CODE RETURNED ',IOS,&
     &   '--> Offending UNIT number: ',LUNIT, &
     &   '--> Name of offending file (if any): ',filename(1:FILELEN), &
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
      !CALL QTRACE(6)
      CALL lsQUIT('ERROR (LSOPEN) TRYING TO OPEN A NON_EXISTING FILE AS OLD',-1)
   else
      write(6,*) 'Wrong error code in error_open! code =', error
   endif

!     List of errors
!
!     Error     Error 
!     Number	Code                            Decription
!     -1	FI_IOSTAT_ENDFIL	        end of file
!     00	FI_IOSTAT_NOTERR	        no error
!     01	FI_IOSTAT_NOTFORSPE	        not FORTRAN-specific error
!     02	FI_IOSTAT_NOTIMP	        not implemented
!     03	FI_IOSTAT_IGNORED	        ignored requested disposition
!     04	FI_IOSTAT_IGNNOTDEL	        ignored requested disposition, file not deleted.
!     17	FI_IOSTAT_SYNERRNAM	        syntax error in NAMELIST input.
!     18	FI_IOSTAT_TOOMANVAL	        too many values for NAMELIST variable
!     19	FI_IOSTAT_INVREFVAR	        invalid reference to variable
!     20	FI_IOSTAT_REWERR	        REWIND error
!     21	FI_IOSTAT_DUPFILSPE	        duplicate file specifications
!     22	FI_IOSTAT_INPRECTOO	        input record too long
!     23	FI_IOSTAT_BACERR	        BACKSPACE error
!     24	FI_IOSTAT_ENDDURREA	        end-of-file during read
!     25	FI_IOSTAT_RECNUMOUT	        record number outside range
!     26	FI_IOSTAT_OPEDEFREQ	        OPEN or DEFINE FILE required
!     27	FI_IOSTAT_TOOMANREC	        too many records in I/O statement
!     28	FI_IOSTAT_CLOERR	        close error
!     29	FI_IOSTAT_FILNOTFOU	        file not found
!     30	FI_IOSTAT_OPEFAI	        open failure
!     31	FI_IOSTAT_MIXFILACC	        mixed file access modes
!     32	FI_IOSTAT_INVLOGUNI	        invalid logical unit number
!     33	FI_IOSTAT_ENDFILERR	        ENDFILE error
!     34	FI_IOSTAT_UNIALROPE	        unit already open
!     35	FI_IOSTAT_SEGRECFOR	        segmented record format error
!     36	FI_IOSTAT_ATTACCNON	        attempt to access non-existent record
!     37	FI_IOSTAT_INCRECLEN	        inconsistent record length
!     38	FI_IOSTAT_ERRDURWRI	        error during write
!     39	FI_IOSTAT_ERRDURREA	        error during read
!     40	FI_IOSTAT_RECIO_OPE	        recursive I/O operation
!     41	FI_IOSTAT_INSVIRMEM	        insufficient virtual memory
!     42	FI_IOSTAT_NO_SUCDEV	        no such device
!     43	FI_IOSTAT_FILNAMSPE	        file name specification error
!     44	FI_IOSTAT_INCRECTYP	        inconsistent record type
!     45	FI_IOSTAT_KEYVALERR	        keyword value error
!     46	FI_IOSTAT_INCOPECLO	        inconsistent OPEN/CLOSE parameters
!     47	FIO_DEF(FI_IOSTAT_WRIREAFIL	write to READONLY file
!     48	FIO_DEF(FI_IOSTAT_INVARGFOR	invalid argument for IO library
!     49	FIO_DEF(FI_IOSTAT_INVKEYSPE	invalid key specification
!     50	FIO_DEF(FI_IOSTAT_INCKEYCHG	inconsistent or duplicate key
!     51	FIO_DEF(FI_IOSTAT_INCFILORG	inconsistent file 
!     52	FIO_DEF(FI_IOSTAT_SPERECLOC	specified record locked
!     53	FIO_DEF(FI_IOSTAT_NO_CURREC	no current record
!     54	FIO_DEF(FI_IOSTAT_REWRITERR	REWRITE error
!     55	FIO_DEF(FI_IOSTAT_DELERR	DELETE error
!     56	FIO_DEF(FI_IOSTAT_UNLERR	UNLOCK error
!     57	FIO_DEF(FI_IOSTAT_FINERR	FIND error
!     59	FIO_DEF(FI_IOSTAT_LISIO_SYN	list-directed I/O syntax error
!     60	FIO_DEF(FI_IOSTAT_INFFORLOO	infinite format loop
!     61	FIO_DEF(FI_IOSTAT_FORVARMIS	format/variable-type mismatch
!     62	FIO_DEF(FI_IOSTAT_SYNERRFOR	syntax error in format
!     63	FIO_DEF(FI_IOSTAT_OUTCONERR	output conversion error
!     64	FIO_DEF(FI_IOSTAT_INPCONERR	input conversion error
!     66	FIO_DEF(FI_IOSTAT_OUTSTAOVE	output statement overflows record
!     67	FIO_DEF(FI_IOSTAT_INPSTAREQ	input statement requires too much data
!     68	FIO_DEF(FI_IOSTAT_VFEVALERR	variable format expression value error
!     70	FIO_DEF(FI_IOSTAT_INTOVF	integer overflow
!     71	FI_IOSTAT_INTDIV	        integer divide by zero
!     72	FI_IOSTAT_FLTOVF	        floating overflow
!     73	FI_IOSTAT_FLTDIV	        floating/decimal divide by zero
!     74	FI_IOSTAT_FLTUND	        floating underflow
!     77	FI_IOSTAT_SUBRNG	        subscript out of range
!     80	FI_IOSTAT_WRONUMARG	        wrong number of arguments
!     81	FI_IOSTAT_INVARGMAT	        invalid argument to math library
!     82	FI_IOSTAT_UNDEXP	        undefined exponentiation
!     83	FI_IOSTAT_LOGZERNEGF	        logarithm of zero or negative value
!     84	FI_IOSTAT_SQUROONEG	        SQUROONEG square root of negative value
!     87	FI_IOSTAT_SIGLOSMAT	        significance lost in math library
!     88	FI_IOSTAT_FLOOVEMAT	        floating overflow in math library
!     89	FI_IOSTAT_FLOUNDMAT	        floating underflow in math library
!     93	FI_IOSTAT_ADJARRDIM	        adjustable array dimension error
!     94	FI_IOSTAT_NEGVEC	        negative vector length in array math function
!     95	FI_IOSTAT_DOMERR	        invalid argument to array math function
!     96	FI_IOSTAT_OVEEXE	        floating point overflow in array math function
!     97	FI_IOSTAT_SIGLOS	        loss of significance in array math function
!     98	FI_IOSTAT_DENNUM	        denormalized floating pointnumber detected in array math function
!     99	FI_IOSTAT_NOTCM   	        Formatted IO is not supported on PFS files
!
end subroutine error_open

!> \brief Routine for opening files in the DALTON program.
!> \author S. Host
!> \date June 2010
!>
!> Based on LSCLOSE by K. Ruud
!> Purpose:
!>   General purpose routine for closing files in the Dalton program.
!>   The routine will ensure that files that are closed and that will not
!>   be reopened in a later part of the program will make their unit 
!>   number available for reuse.
!>
!>   These files are strongly machine dependent, although care has been 
!>   taken to avoid using unit numbers that are illegal or reserved on a
!>   particular architecture. However, the use of LSOPEN and LSCLOSE should
!>   remove much of the machine dependence in the rest of the Dalton 
!>   program
!>
!> Input:
!>  LUNIT    Mandatory
!>  STATUS   Indicate whether the file should be removed or kept
!>           ('KEEP' or 'DELETE'). Note that the default is that the file
!>           is to be deleted
!>
subroutine lsclose(LUNIT,DISP)
implicit none
     !> logical unit number of file to be closed
     integer, intent(inout) :: lunit
     !> Status - 'KEEP' or 'DELETE'
     character(len=*) :: disp
     logical :: fileexists, fileopen,file_exsist
     character(len=200) :: returnfilename
     integer :: ios
!
!     We first deal with the unit number
!
      IF ((LUNIT < 1) .or. (LUNIT > 99) .or. &
     &    (LUNIT == 5) .or. (LUNIT == 6)) call error_close(9001,lunit)
!
!     Check that the file actually has been opened
!
      IF (IUNTAB(LUNIT) == 0) call error_close(9002,lunit)

      INQUIRE(UNIT=LUNIT,EXIST=fileexists,OPENED=fileopen)

      IF (.not. fileexists .or. .not. fileopen) THEN
         call lsquit('LSCLOSE: Tried to close a non-existent or already closed file',-1)
     !    WRITE (LUPRI,'(/A/A,I3)')
     !&     ' WARNING: Tried to close a non-existent or already '//
     !&     'closed file', ' Unit number was: ',LUNIT
     !    CALL FLSHFO(LUPRI)
      END IF

      IF (DISP == 'KEEP') THEN
         CLOSE (LUNIT,STATUS='KEEP',iostat=ios)
         if (ios /= 0) call error_close(9003,lunit)

         IUNTAB(LUNIT) = 0
         LUNIT = -10000-LUNIT
      else if (disp == 'DELETE') THEN
         CLOSE (LUNIT,STATUS='DELETE',iostat=ios)
         if (ios /= 0) then
            INQUIRE(LUNIT,EXIST=file_exsist,NAME=returnfilename,OPENED=fileopen)
            if(ios.GT.0)then
               print*,'ios              = ',ios
               print*,'file_exsist      = ',file_exsist
               print*,'NAME of filename = ',returnfilename
               print*,'file opened      = ',fileopen
            else
               print*,'end-of-file condition is encountered and no error condition is encountered'
            endif
               call lsquit('Something wrong when closing/deleting file',-1)
         endif
         IUNTAB(LUNIT) = 0
         LUNIT = -20000-LUNIT
      else
         call lsquit('Unknown DISP in LSCLOSE',-1)
      END IF
!
!     We release that saved unit number by resetting it to -(code)-LUNIT
!

end subroutine lsclose

subroutine error_close(error,lunit)
implicit none
     integer, intent(in) :: error, lunit

 if (error == 9001) then
      WRITE (6,'(//A/A,I15/A//A)') &
     &   '--> ERROR (LSCLOSE) TRYING TO CLOSE AN ILLEGAL FILE NUMBER', &
     &   '--> SOMEBODY HAS TRIED TO CLOSE UNIT',LUNIT, &
     &   '--> THE PROGRAM DOES NOT ALLOW THE USE OF THIS UNIT NUMBER', &
     &   '### Please report the problem to dalton-admin@kjemi.uio.no' 
      !CALL QTRACE(6)
      CALL lsQUIT('ERROR (LSCLOSE) ILLEGAL FILE NUMBER REQUESTED',-1)
 else if (error == 9002) then
      WRITE (6,'(//A/A,I15//A)') &
     &   '--> ERROR (LSCLOSE) TRYING TO CLOSE A FILE NOT IN USE', &
     &   '--> SOMEBODY IS TRYING TO USE A FILENUMBER THAT HAS '// &
     &   'NOT BEEN USED :',LUNIT, &
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
      !CALL QTRACE(6)
      CALL lsQUIT('ERROR (LSCLOSE) TRYING TO CLOSE A FILE NOT IN USE',-1)
 else if (error == 9003) then
      WRITE (6,'(//A/A/A//A)') &
     &   '--> ERROR (LSCLOSE) TRYING TO KEEP A SCRATCH FILE', &
     &   '--> A FILE HAS BEEN INDICATED TO BE CLOSED AND KEPT,', &
     &   '--> BUT IT APPEARS THE FILE IS A SCRATCH FILE', &
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
      !CALL QTRACE(6)
      CALL lsQUIT('ERROR (LSCLOSE) TRYING TO CLOSE A SCRATCH FILE ',-1)
 else
    call lsquit('Wrong error code in error_close',-1)
 endif
end subroutine error_close

end module files
