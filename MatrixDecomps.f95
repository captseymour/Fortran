module MatrixDecomps

implicit none

contains

!****************************************************************

  subroutine ludcmp(A,L,U,P)

! This subroutine does an LU decomposition using Crout's algorithm,
! with (partial) pivoting.  The matrix is assumed to be square, nxn

    implicit none
    real, intent(in), dimension(:,:) :: A
    real, intent(out), dimension(:,:) :: L,U,P
    real, dimension(:,:), allocatable :: LU
    real, dimension(:), allocatable :: swp
    integer :: pivot(1)
    integer :: i,j,n,k,kp
    real :: csum,rsum

! Initialize.  I can't change A - it may be needed elsewhere in
! the main program, so I'll just load it into a workspace that I
! can change here.  Pivoting will be easier to think about this way.

    n=size(A,1)
    allocate(LU(n,n),swp(n))
    LU=A
    P=0.0
    do i=1,n
      P(i,i)=1.0
    end do

! Is there a need to pivot?  The array output of MAXLOC identifies 
! ROW, COL.  When we pivot, we need to keep track of the pivoting
! in the permutation matrix P.

    pivot=maxloc(LU(:,1))
    if (pivot(1) /= 1) then
      swp=LU(pivot(1),:)
      LU(pivot(1),:)=LU(1,:)
      LU(1,:)=swp
      swp=P(pivot(1),:)
      P(pivot(1),:)=P(1,:)
      P(1,:)=swp
    end if

! Then the first row is unchanged and becomes U's.  The first
! column needs to be finished to become L's.

    LU(2:n,1)=LU(2:n,1)/LU(1,1)

! Now move into the body of the workspace.

    do k=2,n-1

! From the kth pivot element, do the column arithmetic without 
! dividing just yet.    

      do i=k,n
        csum=dot_product(LU(i,1:k-1), LU(1:k-1,k))
        LU(i,k)=LU(i,k)-csum
      end do

! Do we need to pivot?

      pivot=maxloc(LU(k:n,k))
      if (pivot(1) /= k) then
        swp=LU(pivot(1),:)
        LU(pivot(1),:)=LU(k,:)
        LU(k,:)=swp
        swp=P(pivot(1),:)
        P(pivot(1),:)=P(k,:)
        P(k,:)=swp
      end if

! Finish the columns with division.

      LU(k+1:n,k)=LU(k+1:n,k)/LU(k,k)

! Finish the rows with the arithmetic.

      do j=k,n
        rsum=dot_product(LU(k,1:k-1), LU(1:k-1,j))
        LU(k,j)=LU(k,j)-rsum
      end do

    end do

! Take care of that last element of U, where k=n.

    LU(n,n)=LU(n,n)-dot_product(LU(n,1:n-1), LU(1:n-1,n))

! Finally, sort out what is L and what is U.

    L=0.0
    do i=1,n
      L(i,i)=1.0
      L(i+1:n,i)=LU(i+1:n,i)
    enddo
    U=0.0
    do j=1,n
      U(j,j:n)=LU(j,j:n)
    end do

! Clean up before leaving...

    deallocate(LU,swp)

  end subroutine ludcmp

!****************************************************************

  subroutine BackSub(U,v,g)

! This subroutine does backward substitution involving the upper-
! triangular matrix U and the vector v, and returns the result in
! the vector g.

    implicit none
    real, intent(in), dimension(:,:) :: U
    real, intent(in), dimension(:) :: v
    real, intent(out), dimension(:) :: g
    real, allocatable, dimension(:) :: dumv
    integer :: i,n

    n=size(v)
    allocate(dumv(n))
    g=v*0.0
    dumv=v

    do i=n,2,-1
      g(i)=dumv(i)/U(i,i)
      dumv(1:(i-1)) = dumv(1:(i-1))-U(1:(i-1),i)*g(i)
    enddo
    g(1)=dumv(1)/U(1,1)
    
    deallocate(dumv)

  end subroutine BackSub

!****************************************************************

  subroutine ForwSub(L,v,g)

! This subroutine does forward substitution involving the lower-
! triangular matrix L and the vector v, and returns the result in
! the vector g.

    implicit none
    real, intent(in), dimension(:,:) :: L
    real, intent(in), dimension(:) :: v
    real, intent(out), dimension(:) :: g
    real, allocatable, dimension(:) :: dumv
    integer :: i,n

    n=size(v)
    allocate(dumv(n))
    g=v*0.0
    dumv=v

    do i=1,n-1
      g(i)=dumv(i)/L(i,i)
      dumv((i+1):n) = dumv((i+1):n)-L((i+1):n,i)*g(i)
    enddo
    g(n)=dumv(n)/L(n,n)
    
    deallocate(dumv)

  end subroutine ForwSub

end module MatrixDecomps
