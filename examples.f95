program examples

! This program, in addition to being an example of a simple
! Fortran program, specifically illustrates the use of functions
! and subroutines.
!
! The function in this program will sum 1/(n**2) for the first 
! 1000 values of n.
!
! The subroutine will evaluate the quadratic formula for an 
! an input set of coefficients.

implicit none
real :: answer                        !The function result
real, dimension(3) :: coeffs          !The quadratic coefficients.
complex :: root1,root2                !The roots of the quadratic
real :: pi
integer :: i

! First, the function.

answer = sumn2inv(1000)
print*, answer," is from the function"
pi=4.0*atan(1.0)                      !This is a quick way to get pi.
print*, "It should be close to ",pi*pi/6.0

! Next, the subroutine.

print*, "Input the coefficients of a quadratic equation:"
read*, (coeffs(i),i=1,3)

call roots(coeffs,root1,root2)
print*, "The roots are ",root1," and ",root2,"."

contains

!*********************************************************************

  function sumn2inv(nmax) result(sumit)

! This function sums 1/n**2 from 1 to nmax.  NOTE:  SUM is an intrinsic
! function, so you can't use it as a variable name!
!
! Input:  nmax
! Output:  sumit
! Mnemonic:  SUMs N**2 INVerse

    integer, intent(in) :: nmax          !The input upper bound.
    real :: sumit                        !The function variable.
    real, dimension(:), allocatable :: ninv
                                         !Holds 1/n, n=1,nmax
    integer :: n

    allocate(ninv(nmax))                 !Make space for ninv.

    ninv = (/ (1.0/real(n), n=1,nmax) /)         !Create ninv.
    sumit = dot_product(ninv,ninv)       !Sum the squares of ninv.

    deallocate(ninv)                     !Reclaim the space.

  end function sumn2inv

!*********************************************************************

  subroutine roots(a,r1,r2)

! This subroutine returns the roots r1 and r2 of a quadratic equation
! a1*x**2 + a2*x + a3.
!
! Input:  3-dimensional vector a
! Output:  two complex roots, r1 and r2
! Mnemonic: ROOTS of a quadratic equation

    real, dimension(3), intent(in) :: a  !Quadratic coefficients.
    complex, intent(out) :: r1,r2        !Roots, possibly complex.
    complex :: dummy                     !Under the radical.

    dummy = a(2)*a(2) - 4.0*a(1)*a(3)

! I need to use a complex dummy variable above because the square
! root of a real is a real, and the square root of a complex is a
! complex.  The square root of a negative real number, then, gives
! an error; the square root of a negative complex number does not.

    r1 = ( -a(2) + sqrt(dummy) )/2.0/a(1)
    r2 = ( -a(2) - sqrt(dummy) )/2.0/a(1)

  end subroutine roots

end program examples

