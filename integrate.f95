!***********************************************************************
! Integration Module for STAT 8060 Assignment #2.
!***********************************************************************

module integrate

! This module contains routines that perform numerical integration.
! With the exception of importance sampling, these routines are meant 
! to be as general as possible.

  implicit none

contains

!***********************************************************************

  subroutine importance(f)

! This subroutine generates the data for Question 9, which has 
! you study the effect of the importance function when you're
! integrating via Importance Sampling.

    integer, parameter :: sampmax=1000,simmax=5000
    interface
      function f(x)
        real :: f
        real, intent(in) :: x
      end function
    end interface
    integer :: iseed(12),i,j
    real, dimension(sampmax) :: sample,summand
    real :: pi,gxj,integral

! Setup.

    iseed=630507
    call random_seed(put=iseed)
    pi=4.0*atan(1.0)

! Open and designate the output files.

    open(unit=10,file='isuni.dat')
    open(unit=11,file='iscau.dat')
    open(unit=12,file='ist30.dat')

! Generate 5000 data points in each instance.

    do i=1,simmax

! For the U(-50,50) importance function...

      call random_number(sample)
      sample=sample*50.0
      summand=0.0
      do j=1,sampmax

! The U(-50,50) density is g(x) = 1/100...

        summand(j)=f(sample(j))*100.0
      end do
      integral=sum(summand)/sampmax
      write(10,*) integral

! For the Standard Cauchy importance function...

      call random_number(sample)
      sample=tan(pi*(sample-0.5))
      do
        where (abs(sample) > 50.0)
          sample = tan(pi*(sample-0.5))
        end where
        if (maxval(abs(sample)) < 50.0) exit
      end do

      summand=0.0
      do j=1,sampmax

! Because the Standard Cauchy density is g(x) = 1/(pi*(1+x*x))...

        gxj=1.0/pi/(1.0+sample(j)*sample(j))
        summand(j)=f(sample(j))/gxj
      end do
      integral=sum(summand)/sampmax
      write(11,*) integral

! Now for the t(30) importance function...

      sample=tdist(30,sampmax)
      do
        where (abs(sample) > 50.0)
          sample = tdist(30,1)
        end where
        if (maxval(abs(sample)) < 50.0) exit
      end do

      do j=1,sampmax

! Here's the t density.  

        gxj=gamma(15.5)/sqrt(pi*30.0)/gamma(15.0)
        gxj=gxj/(1.0+sample(j)*sample(j)/30)**15.5
        summand(j)=f(sample(j))/gxj
      end do
      integral=sum(summand)/sampmax
      write(12,*) integral

    end do

    close(10)
    close(11)
    close(12)

    contains

      function tdist(df,n) result(tvec)

! I originally wrote this function to return an n-dimensional 
! vector, but the code to generate the sample runs faster if 
! this returns one scalar at a time.

        integer, intent(in) :: df,n
        real, dimension(n) :: tvec
        real, dimension(n*(df+1)) :: unif,normal
        integer :: lgth,i,j
        real :: pi,den

        pi=4.0*atan(1.0)
        lgth=n*(df+1)
        
! Generate the normals we'll need to build the t.

        call random_number(unif)
        normal(1:lgth/2) = sqrt( -2.0*log(unif(1:lgth/2)) )* &
                           cos( 2.0*pi*unif(lgth/2+1:lgth) )
        normal(lgth/2+1:lgth) = sqrt( -2.0*log(unif(1:lgth/2)) )* &
                                sin( 2.0*pi*unif(lgth/2+1:lgth) )

! Initialize the loop.

        i=1
        j=1
	tvec=0.0
        do

! Make the Chi-Square variate for the denominator.
          
          den=dot_product(normal(i:i+df),normal(i:i+df))

! Now use the final normal to create the t variate.

          tvec(j)=normal(i+df+1)/sqrt(den/df)
          
! Update the counters and quit if we're done.

          i=i+df
          j=j+1
          if (j > n) exit
        end do

      end function tdist

  end subroutine importance

!***********************************************************************

  function midpoint(f,a,b,tol)

! This function uses the extended midpoint rule to evaluate
! the integral of the function given in F.  The only input to
! this function is a value X.

    integer, parameter :: nmax=10
    real :: midpoint
    interface
      function f(x)
        real :: f
        real, intent(in) :: x
      end function
    end interface
    real, intent(in) :: a,b,tol
    integer :: n
    real :: oldI,newI

    print*, "Midpoint Rule Iterates"

! Initialize

    oldI = midn(f,1,a,b)
    print'(i2,1x,f12.9)', 1,oldI

! Iterate until convergence or until we're sick of it.

    do n=2,nmax

! Evaluate this iteration.

      newI = midn(f,n,a,b)
      print'(i2,1x,f12.9)', n,newI
  
! Check for convergence

      if (abs(newI-oldI) <= tol*abs(oldI)) then
        midpoint = newI
        return
      endif
  
! If we haven't converged, then let's go around again.

      oldI=newI

    enddo

! If we make it to this step, then we're officially sick of
! iterating.

    midpoint = newI

  end function midpoint

!***********************************************************************

  function midn(f,n,a,b)

! This function evaluates the nth step of the midpoint rule.

    real :: midn
    interface
      function f(x)
        real :: f
        real, intent(in) :: x
      end function
    end interface
    integer, intent(in) :: n
    real, intent(in) :: a,b
    integer :: i
    real :: h,old,midsum,spot
    save old

    h=(b-a)/(3**(n-1))
    if (n == 1) then
      midn = (b-a)*f((a+b)/2.0)
      old=midn
      return
    endif

    midsum=0.0
    do i=1,2*(3**(n-2))-1,2
      spot = 3.0*i*h/2.0
      midsum = midsum + f(a+spot-h) + f(a+spot+h)
    enddo

    midn = old/3.0 + h*midsum
    old=midn

  end function midn

!***********************************************************************

  function mc2anti(f,a,b)

! This function uses Monte Carlo averaging along with 
! antithetic variables to evaluate the integral of the function 
! given in F.  The only input to this function is a value X.

    integer, parameter :: nmax=500000
    real :: mc2anti
    interface
      function f(x)
        real :: f
        real, intent(in) :: x
      end function
    end interface
    real, intent(in) :: a,b
    integer :: idum(12),n
    real, dimension(nmax) :: x,fx

    print*, "Monte Carlo answer after ",nmax," iterates"

! Initialize.

    idum=630507
    call random_seed(put=idum)

! Grab half the U(0,1)'s you need, then complete the antithetic
! sample and transform it all into U(a,b).

    call random_number(x(1:nmax/2))
    x(nmax/2+1:nmax)=1.0-x(1:nmax/2)
    x=x*(b-a)+a
 
! Get f(x).

    do n=1,nmax
      fx(n)=f(x(n))
    enddo

! And compute the estimate.

    mc2anti = (b-a)*sum(fx)/nmax
    print*, mc2anti

  end function mc2anti

!***********************************************************************

  function mc2(f,a,b)

! This function uses Monte Carlo averaging to evaluate the integral 
! of the function given in F.  The only input to this function is a 
! value X.

    integer, parameter :: nmax=500000
    real :: mc2
    interface
      function f(x)
        real :: f
        real, intent(in) :: x
      end function
    end interface
    real, intent(in) :: a,b
    integer :: idum(12),n
    real, dimension(nmax) :: x,fx

    print*, "Monte Carlo answer after ",nmax," iterates"

! Initialize

    idum=650507
    call random_seed(put=idum)

! Grab some U(0,1)'s and transform them to U(a,b)'s.

    call random_number(x)
    x=x*(b-a)+a

! Get f(x)

    do n=1,nmax
      fx(n)=f(x(n))
    enddo

! And compute the estimate

    mc2 = (b-a)*sum(fx)/nmax
    print*, mc2

  end function mc2

!***********************************************************************

  function mc1(f,a,b,bnd)

! This function uses Monte Carlo acceptance-rejection integration 
! to evaluate the integral of the function given in F.  The only 
! input to this function is a value X.

    integer, parameter :: nmax=1e6
    real :: mc1
    interface
      function f(x)
        real :: f
        real, intent(in) :: x
      end function
    end interface
    real, intent(in) :: a,b,bnd
    integer :: idum(12),n
    real, dimension(nmax) :: x,y
    real :: kount,area

    print*, "Monte Carlo Answer after ",nmax," iterations"

! Initialize

    kount=0.0
    area=bnd*(b-a)
    idum=630507
    call random_seed(put=idum)

! Draw a uniform points in the box.

    call random_number(x)
    x=x*(b-a)+a
    call random_number(y)
    y=y*bnd
  
! If it is under the curve, accumulate it.  Otherwise, don't.

    do n=1,nmax
      if (y(n) < f(x(n))) kount=kount+1.0
    end do

! Compute the estimate

    mc1=area*kount/nmax
    print*, mc1

  end function mc1

!***********************************************************************

  function glegquad(f)

! This function uses 16 Gauss-Legendre weights and nodes to compute
! the integral between 0 and 1.  Note that f(x) must be defined
! in such a way that the limits of integration match up with those
! of the chosen Gaussian quadrature.

    real :: glegquad
    interface
      function f(x)
        real :: f
        real, intent(in) :: x
      end function
    end interface
    real, dimension(16) :: x,w,fx
    integer :: i

    x=(/0.0483076656877383162348126, 0.1444719615827964934851864, &
        0.2392873622521370745446032, 0.3318686022821276497799168, &
        0.4213512761306353453641194, 0.5068999089322293900237475, &
        0.5877157572407623290407455, 0.6630442669302152009751152, &
        0.7321821187402896803874267, 0.7944837959679424069630973, &
        0.8493676137325699701336930, 0.8963211557660521239653072, &
        0.9349060759377396891709191, 0.9647622555875064307738119, &
        0.9856115115452683354001750, 0.9972638618494815635449811/)

    w=(/0.0965400885147278005667648, 0.0956387200792748594190820, &
        0.0938443990808045656391802, 0.0911738786957638847128686, &
        0.0876520930044038111427715, 0.0833119242269467552221991, &
        0.0781938957870703064717409, 0.0723457941088485062253994, &
        0.0658222227763618468376501, 0.0586840934785355471452836, &         
        0.0509980592623761761961632, 0.0428358980222266806568786, &
        0.0342738629130214331026877, 0.0253920653092620594557526, &
        0.0162743947309056706051706, 0.0070186100094700966004071/)

! Find f(x) for the nodes x.  I have to loop here because of 
! the way f(.) is defined.

    do i=1,16
      fx(i) = f(x(i))
    end do

    glegquad = dot_product(fx,w)

  end function glegquad

!***********************************************************************

  function romberg(f,a,b,tol)

! This function uses Romberg's integration rule to evaluate
! the integral of the function given in F.  The only input to
! this function is a value X.

    integer, parameter :: nmax=20
    real :: romberg
    interface
      function f(x)
        real :: f
        real, intent(in) :: x
      end function
    end interface
    real, intent(in) :: a,b,tol
    integer :: n,j
    real, dimension(nmax) :: oldI,newI

    write(9,'(a)') "Romberg Iterates"

! Initialize

    oldI(1) = trapn(f,1,a,b)
    print'(i2,1x,f12.9)', 1,oldI(1)

! Iterate until convergence or until we're sick of it.

    do n=2,nmax

! Evaluate this row.

      newI(1) = trapn(f,n,a,b)
      do j=2,n
        newI(j)=(4**(j-1)*newI(j-1) - oldI(j-1))/(4**(j-1)-1)
      enddo  
      print'(i2,10(1x,f12.9))', n,newI(1:n)
  
! Check for convergence

      if (abs(newI(n)-newI(n-1)) <= tol*abs(newI(n-1))) then
        romberg = newI(n)
        return
      endif
  
! If we haven't converged, then let's go around again.

      oldI=newI

    enddo

! If we make it to this step, then we're officially sick of
! iterating.

    romberg = newI(nmax)

  end function romberg


!***********************************************************************

  function trapezoid(f,a,b,tol)

! This function uses the extended trapezoidal rule to evaluate
! the integral of the function given in F.  The only input to
! this function is a value X.

    integer, parameter :: nmax=20
    real :: trapezoid
    interface
      function f(x)
        real :: f
        real, intent(in) :: x
      end function
    end interface
    real, intent(in) :: a,b,tol
    integer :: n
    real :: oldI,newI

    print*, "Trapezoidal Rule Iterates"

! Initialize

    oldI = trapn(f,1,a,b)
    print'(i2,1x,f12.9)', 1,oldI

! Iterate until convergence or until we're sick of it.

    do n=2,nmax

! Evaluate this iteration.

      newI = trapn(f,n,a,b)
      print'(i2,1x,f12.9)', n,newI
  
! Check for convergence

      if (abs(newI-oldI) <= tol*abs(oldI)) then
        trapezoid = newI
        return
      endif
  
! If we haven't converged, then let's go around again.

      oldI=newI

    enddo

! If we make it to this step, then we're officially sick of
! iterating.

    trapezoid = newI

  end function trapezoid

!***********************************************************************

  function trapn(f,n,a,b)

! This function evaluates the nth step of the trapezoidal rule.

    implicit none
    real :: trapn
    interface
      function f(x)
        real :: f
        real, intent(in) :: x
      end function
    end interface
    integer, intent(in) :: n
    real, intent(in) :: a,b
    integer :: i
    real :: old,trapsum,x
    save old

    if (n == 1) then
      trapn = (b-a)*(f(b)+f(a))/2.0
      old=trapn
      return
    endif

    trapsum=0.0
    do i=1,2**(n-2)
      x = a + (2*i-1)*(b-a)/(2**(n-1))
      trapsum = trapsum + f(x)
    enddo

    trapn = 0.5*old + (b-a)*trapsum/(2**(n-1))
    old=trapn

  end function trapn

end module integrate
