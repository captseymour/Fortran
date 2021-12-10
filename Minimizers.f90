module Minimizers

! This module contains routines for the numerical minimization of
! functions.  These routines are meant to be as general as possible.

implicit none

contains

!****************************************************************

subroutine evolve(xmax,func,fit,ngen,maxgen)

! This subroutine uses the evolutionary algorithm to optimize the 
! function given in FUNC.

use random

implicit none
integer, parameter :: nmf=10   !Number of males = number of females.
real, parameter :: eps=1e-6
integer, intent(in) :: maxgen
integer, intent(out) :: ngen
real, intent(out) :: xmax,fit
real, external :: func
integer :: iseed=-12345
integer, dimension(1) :: locmax
real, dimension(nmf,2) :: male,female  !The second component for each
real, dimension (:,:) :: offspring     !is the fitness.
integer :: i,j,index
real :: pi,unif(1),xmax0,fit0
allocatable offspring

! Initialize the subroutine, including the random number generator.

allocate(offspring(nmf*nmf,2))
pi=4.0*atan(1.0)
male=0.0
female=0.0
offspring=0.0
call ran_seed(sequence=iseed)

! Let's watch the population evolve!

open(unit=10,file="evolve.out")

! Initial population.  Start 'em out Cauchy centered about 0.
! The fitness of the initial population is irrelevant in this
! scheme.

write(10,*) "Initial population, with fitness"
do i=1,nmf
  call ran(unif)
  male(i,1)=tan(pi*(unif(1)-0.5))
  call ran(unif)
  female(i,1)=tan(pi*(unif(1)-0.5))
  write(10,*) male(i,1),func(male(i,1)),female(i,1),func(female(i,1))
enddo
write(10,*)

xmax0=-1e9
fit0=-1e9

do ngen=1,maxgen

! Now they get to "do their thing".  Each male gets to reproduce 
! with each female.  The mode of reproduction is a simple average.

  do i=1,nmf
    do j=1,nmf
      index=(i-1)*nmf+j
      offspring(index,1)=(male(i,1)+female(j,1))/2.0

! There is a 10% chance of a mutation.  Does this one get to mutate?

      call ran(unif)
      if (unif(1) < 0.1) then
        call ran(unif)
        offspring(index,1) = offspring(index,1) + (2.0*unif(1)-1.0)
      endif

! And evaluate!

      offspring(index,2)=func(offspring(index,1))
    enddo
  enddo

! Okay, who gets to stay for the next generation?  We'll alternately
! assign males and females, "Ladies First".  :-)

  do i=1,nmf
    locmax=maxloc(offspring(:,2))
    female(i,:)=offspring(locmax(1),:)
    xmax=female(i,1)
    fit=female(i,2)
    offspring(locmax(1),2)=-1e9
    locmax=maxloc(offspring(:,2))
    male(i,:)=offspring(locmax(1),:)
    offspring(locmax(1),2)=-1e9
  enddo

! Immigrate.  The immigrant is Cauchy, centered about the most fit.
! It will replace the least fit - which is the last male chosen.

!  call ran(unif)
!  male(nmf,1)=tan(pi*(unif(1)-0.5))+female(1,1)
!  male(nmf,2)=func(male(nmf,1))

! Check for convergence and go for another generation.

  write(10,*) "Generation ",ngen
  do i=1,nmf
    write(10,*) (male(i,j),j=1,2),(female(i,j),j=1,2)
  enddo
  write(10,*)

  if (abs(fit-fit0) < eps) then
    close(10)
    deallocate(offspring)
    return
  else
    xmax0=xmax
    fit0=fit
  endif

enddo
close(10)
deallocate(offspring)

end subroutine evolve

!************************************************************

subroutine qnmin(n,b,func,f0,h,nevals)

! This subroutine is a quasi-newton minimizer from Nash, "Compact 
! Numerical Methods for Computers" (1979), alg 21.
!
! Modifications:
! B.D. Ripley, 5/1980
! R.L. Smith, 4/1988
! P.L. Seymour 10/1991, 8/1995;
!              changed to Fortran 90 3/2000;
!              reduced to single precision for STAT 8060 11/2001
!
! N   Dimension of parameter B.       
! B   Parameter to be returned which minimizes the function in LKHD.
! P0  Value of pseudo-likelihood at B.
! NEVALS  Maximum number of function evaluations.
!
! Mnemonic:  Quasi-Newton MINimizer.
!
! Subroutines called:  FUNC, GRAD, VECMAT

implicit none

integer, intent(in) :: n,nevals
real, dimension(:), intent(inout) :: b
real, intent(out) :: f0
real, dimension(n,n), intent(out) :: h
real, dimension(:), allocatable :: b0,g0,g,t
real, dimension(:,:), allocatable :: vm1,vm2,vm3
real, external :: func
integer :: count,ifn,i
real :: k,sn,d1,f,d2
real, parameter :: w=0.2,tol=1.e-6,eps=1.e-6

allocate(b0(n),g0(n),g(n),t(n),vm1(n,n),vm2(n,n),vm3(n,n))

! Is this an infeasible point?  If so, exit.

f0=func(n,b)
if (f0 >= 1.e9) then
  stop "QNMIN: Initial point infeasible"
endif

! Otherwise, calculate the gradient and increment IFN, the number of 
! function evaluations, and IG, the number of gradient calculations.

call grad(n,b,f0,func,g)
ifn=n+1

do

! Reset inverse-Hessian.
  h=0.0
  do i=1,n
    h(i,i)=1.0
  enddo

! Top of iteration.
  do
  
! Store parameter and gradient.
    b0=b
    g0=g

! find search direction t
    t=-matmul(h,g)
    sn=dot_product(t,t)
    d1=-dot_product(t,g)

! Check if downhill.
    if(d1 <= 0.0) then
      deallocate(b0,g0,g,t,vm1,vm2,vm3)
      return
    endif
    if(d1 == 0.0) exit

! Search along T...
    sn=0.5/sqrt(sn)
    k=min(1.0,sn)

! ...for a good next step.
    do
      count=0
      b=b0+k*t
      do i=1,n
        if(abs(k*t(i)) < eps) count=count+1
      enddo

! Check if converged.  
      f=func(n,b)
      if(count == n .and. f < 1.e9) then

! Get outta here - you've converged, and you're going to crash
! the inverse-Hessian stuff if you compute it with this f and its
! gradient!  The last inverse-Hessian is the best approximation
! you're going to get at this point.

        f0=f
        deallocate(b0,g0,g,t,vm1,vm2,vm3)
        return
      endif

! Otherwise . . .
      ifn=ifn+1
      if(ifn >= nevals) then
        stop "QNMIN: Too many function evaluations"
      endif
      if(f < f0-d1*k*tol) exit
      k=w*k
    enddo

! New lowest value.
    f0=f
    call grad(n,b,f,func,g)
    ifn=ifn+n

! Update inverse-Hessian.  Re-use some variables to save memory.

    t=k*t     !Change in parameter from last time to this.
    g0=g-g0     !Change in gradient from last time to this.
    d1=dot_product(t,g0)

! check if +ve def addition.  If not, reset the inverse-Hessian.
    if(d1 <= 0.0) exit
    b0=matmul(h,g0)
    d2=dot_product(b0,g0)
    d2=1+d2/d1
    call vecmat(n,t,t,vm1)
    call vecmat(n,t,b0,vm2)
    call vecmat(n,b0,t,vm3)
    h=h+(d2*vm1-vm2-vm3)/d1

! Top of iteration.
  enddo
  
! Reset inverse-Hessian.
enddo

end subroutine qnmin

!************************************************************

subroutine grad(n,b,f,func,g)

! This subroutine approximates a gradient vector for the function
! computed in the subroutine LKHD.
!
! Mnemonic:  GRADient
!
! Called from:  QNMIN
!
! Subroutines called:  FUNC

implicit none
integer, intent(in) :: n
real, intent(in) :: f
real, dimension(:), intent(in) :: b
real, external :: func
real, dimension(n), intent(out) :: g
real, dimension(:), allocatable :: x
real, parameter :: h=1.e-3
integer :: i
real :: fh

allocate(x(n))
x=b

do i=1,n
  x(i)=x(i)+h
  fh=func(n,x)
  g(i)=(fh-f)/h
  x(i)=b(i)
enddo

deallocate(x)
end subroutine grad

!************************************************************

subroutine vecmat(n,a,b,v)

! This subroutine returns the product of an nx1 vector with
! a 1xn vector to yield a rank-1 matrix.

implicit none
integer, intent(in) :: n
real, dimension(:), intent(in) :: a,b
real, dimension(:,:), intent(out) :: v
integer :: i

do i=1,n
  v(i,1:n)=a(i)*b(1:n)
enddo

end subroutine vecmat

end module Minimizers
