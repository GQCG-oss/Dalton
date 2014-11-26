!> @file 
!> @brief Timing module for PBC-SCF.
!>
!> Generic module to manage timings of diffenent tasks. Add times to the
!> diffenent taskt f.eks. in subsequent iterations.
!>
!> Usage: 
!> type(timing_info), pointer :: timings(:)
!> pbc_timings_init(timings, num_jobs)
!> pbc_timings_name_timer(timings, num_of_timer, 'exchange')
!> ..
!> ..
!> !add time num of timer
!> pbc_timings_start(timings, num_of_timer)
!> !execute code to be timed ...
!> pbc_timings_stop(timings, num_of_timer) 
!> ...
!> pbc_timings_print(timings)
!> pbc_timings_destruct(timings)

module pbc_timings
use precision

type timing_info
		integer :: num_jobs
		real(realk) :: cputime
		real(realk) :: walltime
		character (len=20) :: jobname
		real(realk) :: cputime_begin, walltime_begin
end type

contains
	
!> @brief Initialize and allocate timing_info.
!> @param timings Timing_info pointer.
!> @param num_jobs_inp Number of different jobs to time.
subroutine pbc_timings_init(timings, num_jobs_inp)
   implicit none
	! input
	integer, intent(in) :: num_jobs_inp
	type(timing_info), pointer, intent(inout):: timings(:)
	! local
	integer :: i

	allocate(timings(num_jobs_inp))
	
	do i = 1, num_jobs_inp
		timings(i)%cputime = 0.0_realk
		timings(i)%walltime = 0.0_realk
		timings(i)%jobname = '--not initialized --'
		timings%num_jobs = num_jobs_inp
	end do	

end subroutine pbc_timings_init

!> @brief Deallocate mem.
!> @param timings Timing_info pointer.
subroutine pbc_timings_destruct(timings)
		implicit none
		! input
		type(timing_info), pointer, intent(inout):: timings(:)

		deallocate(timings)

end subroutine pbc_timings_destruct

!> @brief Initialize names of the different jobs.
!> @param timings Timing_info pointer.
!> @param jname Jobname.
!> @param job The number of the job to name.
subroutine pbc_timings_name_timer(timings, jname, job)
		implicit none
		! input
		type(timing_info), pointer, intent(inout):: timings(:)
		integer, intent(in) :: job
		character (len=20), intent(in) :: jname

		timings(job)%jobname = jname

end subroutine pbc_timings_name_timer 

!> @brief Add time to one of the job. Start timing.
!> @param timings Timing_info pointer.
!> @param job Number of the job.
subroutine pbc_timings_start(timings, job)
   implicit none
	! input 
   integer, intent(in) :: job
	type(timing_info), pointer, intent(inout) :: timings(:)
	! local
	real(realk) :: cputime, walltime

   !$OMP CRITICAL (timematop)
   call ls_gettim(cputime,walltime)
	timings(job)%walltime_begin = walltime
	timings(job)%cputime_begin = cputime
   !$OMP END CRITICAL (timematop)

end subroutine pbc_timings_start

!> @brief Add time to one of the job. Stop timing.
!> @param timings Timing_info pointer.
!> @param job Number of the job.
!> @param printlevel_i 0 do not print, 1 print to lupri, 2 print to screen and lupri.
!> @param lupri Logical print unit.
subroutine pbc_timings_stop(timings, job, printlevel_i, lupri)
   implicit none
	! input
   integer, intent(in) :: job
	type(timing_info), pointer, intent(inout) :: timings(:)
	integer, intent(in), optional :: printlevel_i
	integer, intent(in), optional :: lupri
	! local
   real(realk) :: cputime_end, walltime_end, deltawall, deltacpu
	integer :: printlevel
	character (len=9) :: dw_s, dc_s
	
	if ( present(printlevel_i) ) then
		printlevel = printlevel_i
	else
		printlevel = 1
	endif

	if (printlevel > 2) then
		if ( .not. present(lupri) ) then
			printlevel = 1
		endif
		if ( lupri == 6 ) then
			printlevel = 1
		endif
	endif

   call ls_gettim(cputime_end,walltime_end)
	deltawall = walltime_end - timings(job)%walltime_begin
	deltacpu = cputime_end - timings(job)%cputime_begin
	
	timings(job)%cputime = timings(job)%cputime + deltacpu
	timings(job)%walltime = timings(job)%walltime + deltawall

	if ( printlevel > 0 ) then
		write (dw_s,'(2F9.4)') deltawall
		write (dc_s,'(2F9.4)') deltacpu
		write (6, *) '>>>> Time consumption: ', &
 				& timings(job)%jobname, ',  Walltime:', dw_s, ',  CPUtime:', dc_s, '.'
	endif
	if ( printlevel == 2 ) then 
		write (dw_s,'(2F9.4)') deltawall
		write (dc_s,'(2F9.4)') deltacpu
		write (lupri, *) '>>>> Time consumption: ', &
 				& timings(job)%jobname, ',  Walltime:', dw_s, ',  CPUtime:', dc_s, '.'
	endif

end subroutine pbc_timings_stop

	
!> @brief Print all the collected information.
!> @param timings Timing_info pointer.
!> @param lupri Print unit.
subroutine pbc_timings_print(timings, lupri)
	implicit none
	! input
	type(timing_info), pointer, intent(in) :: timings(:)
	integer, intent(in) :: lupri
	!local
	integer :: num_jobs, i

	write (lupri, *)
	write (lupri, *) '=================================================================='
	write (lupri, *) '======================= PBC - SCF TIMINGS ========================'
	write (lupri, *) '=================================================================='
	write (lupri,'(A31,A18,A18)') '   Task                      ', &
		& 'Walltime    ', 'Cputime     '
	write (lupri, *) '------------------------------------------------------------------'
	num_jobs = timings(1)%num_jobs
	do i = 1, num_jobs
		write (lupri,'(A5,A20,2F18.4,2F18.4)') ' >>> ', timings(i)%jobname, &
			& timings(i)%walltime, timings(i)%cputime
	enddo
	write (lupri, *) '=================================================================='
	write (lupri, *) '=================================================================='
	write (lupri, *)
	
end subroutine pbc_timings_print

end module pbc_timings
