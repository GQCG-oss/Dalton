module minimize_module
! original code by Stinne Host
! modifications by B. Jansik, Arhus, Jun 2006
use Matrix_module
#define VERBOSE 0

private::SIMPLEX_solver,HIGH_LOW,VECNORM

contains
!Finding the minimum of an n-dimensional function using the Downhill Simplex method.
!func is name of function to be minimized
subroutine simplex(func, n, initial_simplex, solution, params)
Implicit none

     real(realk), external :: func
     integer, intent(in)   :: n   !Dim. of problem
     real(realk), dimension(n+1,n), intent(in) :: initial_simplex
     integer                       :: M 
     real(realk)                   :: tolf, tolx
     real(realk), dimension(n), intent(out) :: solution
     real(realk), dimension(:) :: params

   !Tolerance on f(x1,x2) and simplex:
   tolf = 0.0000001 ; tolx = 0.00001

   call SIMPLEX_solver(func, initial_simplex, n, tolf, tolx, solution, M, params)

#if VERBOSE
print *, "Total number of iterations:", M
print "(/)"
print '("Tolerance on function value:", e10.2)', tolf
print '("Tolerance on Simplex size:  ", e10.2)', tolx
#endif

end subroutine simplex

subroutine SIMPLEX_solver(func, sim, n, tolf, tolx, x_final, M, params)
     real(realk), external :: func
     integer, intent(in) :: n 
     integer, intent(out) :: M
     integer :: i, j, t, RF, E, C, RD
     real(realk) :: p_m(n-1,n), p_new(n+1,n), diff_f, diff_x
     real(realk), dimension(n) :: p_h, p_l, p0, p_r, p_e, p_c,diffp0ph
     real(realk), intent(in) ::  sim(n+1,n), tolf, tolx
     real(realk), intent(out) :: x_final(n)
     real(realk), dimension(:) :: params

RF = 0   !Counters for number of reflections, expansions...
E = 0
C = 0
RD = 0

M = 0    !Counter for number of iterations.

#if VERBOSE
call MATOUT("Initial simplex:", sim, n+1, n)
#endif

!Find highest point (p_h),lowest point (p_l) and centroid (p0) of simplex:

call HIGH_LOW(func, n, sim, p_h, p_l, p_m, p0, params)
 
!Now rolling downhill towards the minimum:

do
   M = M + 1
   if (M == 1000) then
#if VERBOSE
      print *, "Too many iterations!"
#endif
      STOP
      exit
   endif
   diff_f = ABS(func(p_h, params)-func(p_l, params))                !Convergence criteria.
   diffp0ph = p0-p_h
   call VECNORM(diffp0ph, n, diff_x)
   diff_x = ABS(diff_x)
   if (diff_f < tolf .AND. diff_x < tolx) then
#if VERBOSE
      print *, "Convergence!"
      print "(/)"
#endif
      exit
   endif
   t = 1
   p_r = p0 - (p_h - p0)                   ! Try reflection.
   if (func(p_r,params) < func(p_l,params)) then
      E = E + 1
      p_e = 2*p_r - p0                     !Reflection excellent - Expansion.
      p_new(1,1:n) = p_r
      p_new(2,1:n) = p_l
      do i = 3, n+1
         p_new(i,1:n) = p_m(t,1:n)
         t = t + 1
      enddo
      call HIGH_LOW(func, n, p_new, p_h, p_l, p_m, p0,params)
      cycle
   else
      if (func(p_r,params) < func(p_h,params)) then            !Reflection ok.
         RF = RF + 1
         p_new(1,1:n) = p_r
         p_new(2,1:n) = p_l
         do i = 3, n+1
            p_new(i,1:n) = p_m(t,1:n)
            t = t + 1
         enddo
         call HIGH_LOW(func, n, p_new, p_h, p_l, p_m, p0,params)
         cycle
      else
         p_c = p0 + 0.5*(p_h - p0)     !Reflection no good - Contraction.
         C = C + 1
         if (func(p_c,params) < func(p_h,params)) then
            p_new(1,1:n) = p_c
            p_new(2,1:n) = p_l
            do i = 3, n+1
               p_new(i,1:n) = p_m(t,1:n)
               t = t + 1
            enddo
            call HIGH_LOW(func, n, p_new, p_h, p_l, p_m, p0,params)
            cycle
         else
            RD = RD + 1
            p_h = 0.5*(p_h + p_l)   !Contraction no good - Reduction.
            p0 = 0.5*(p0 + p_l)
            do i = 1, n-1
               p_m(i,1:n) = 0.5*(p_m(i,1:n) + p_l)
            enddo
         endif
      endif
   endif
enddo

#if VERBOSE
print *, "Number of operations used:"
print *, "Reflections: ", RF
print *, "Expansion:   ", E
print *, "Contractions:", C
print *, "Reductions:  ", RD
print "(/)"
#endif 
x_final = p_l

end subroutine SIMPLEX_solver

subroutine HIGH_LOW(func, n, POINTS, p_h, p_l, p_m, p0,params)
implicit none

     real(realk), external :: func
     integer :: t, i
     integer, intent(in) :: n
     real(realk), intent(in) :: POINTS(n+1,n)
     real(realk), intent(out) :: p_h(n), p_l(n), p_m(n-1,n), p0(n) 
     real(realk) :: f_points(n+1)
     real(realk), dimension(:) :: params


do i = 1, n+1
   f_points(i) = func(POINTS(i,1:n),params)
enddo

t = 1
                                                                                
do i = 1, n+1
   if (f_points(i) == MAXVAL(f_points)) then
      p_h = POINTS(i,1:n)
   else
      if (f_points(i) == MINVAL(f_points)) then
         p_l = POINTS(i,1:n)
      else
         p_m(t,1:n) = POINTS(i,1:n)
         t = t +1
      endif
   endif
enddo

do i = 1, n
   p0(i) = (1.0/n)*(SUM(POINTS(1:n+1,i)) - p_h(i))
enddo

end subroutine HIGH_LOW

!Calculation of the norm-square of a vector c of dimension m.
subroutine VECNORM(c, m, N)
  implicit none
  integer,intent(in) :: m
  real(realk),intent(in) :: c(m)
  real(realk),intent(inout) :: N
  real(realk), external :: ddot
  
  N=SQRT(ddot(m,c,1,c,1))
end subroutine VECNORM

end module minimize_module


