! The following module is a user-written Fortran 90 uniform random 
! number generator.  Its use is very similar to the intrinsic Fortran 
! 90 functions RANDOM_SEED and RANDOM_NUMBER.
!
! To use the random number generator in your own program, follow the
! instructions using this MODULE, named RANDOM, given in the Tutorial
! on the course web page.  
!
! You may also use this MODULE within your program, but it will have
! to violate the good programming manners that I'm trying to teach 
! you.  You will have to put this entire MODULE before your PROGRAM 
! statement.  

MODULE random

! Adapted from _Numerical Recipes in Fortran 90_
!
! This module supports the random number routine RAN.  It provides 
! the generator with five vectors of integers for use as internal 
! state space.  The first three integers (IRAN, JRAN, KRAN) are 
! maintained as nonnegative values, while the last two (MRAN, NRAN) 
! have 32-bit nonzero values.  Also provided by this module is support 
! for initializing or reinitializing the state space to a desired 
! standard seque3nce number, hashing the initial values to random 
! values, and allocating and deallocating the internal workspace.

  implicit none

!Kind values for (ideally) 32-bit integers and single precision.
  integer, parameter :: k4b=selected_int_kind(9),sp=kind(1.0)
  integer(kind=k4b), parameter :: hg=huge(1_k4b),hgm=-hg,hgng=hgm-1
  integer(kind=k4b), save :: lenran=0,seq=0
  integer(kind=k4b), save :: iran0,jran0,kran0,nran0,mran0,rans
  integer(kind=k4b), dimension (:,:), pointer, save :: ranseeds
  integer(kind=k4b), dimension (:), pointer, save :: iran,jran,kran,mran,nran,ranv
  real(kind=sp), save :: amm

CONTAINS

  FUNCTION arth(first,increment,n) result(arth_result)
    integer(kind=k4b), intent(in) :: first, increment, n
    integer(kind=k4b), dimension(n) :: arth_result
    integer(kind=k4b) :: k
    if (n>0) arth_result(1)=first
    do k=2,n
      arth_result(k)=arth_result(k-1)+increment
    enddo
  END FUNCTION arth

  FUNCTION reallocate_iv(p,n) result(real_iv_res)
    integer(kind=k4b), dimension (:), pointer :: p, real_iv_res
    integer(kind=k4b), intent(in) :: n
    integer(kind=k4b) :: nold,ierr
    allocate(real_iv_res(n),stat=ierr)
    if(ierr /= 0) stop "REALLOCATE_IV: problem in attempt to allocate memory"
    if (.not. associated(p)) return
    nold=size(p)
    real_iv_res(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  END FUNCTION reallocate_iv

  FUNCTION reallocate_im(p,n,m) result(real_im_res)
    integer(kind=k4b), dimension(:,:), pointer :: p, real_im_res
    integer(kind=k4b), intent(in) :: n,m
    integer(kind=k4b) :: nold, mold, ierr
    allocate(real_im_res(n,m),stat=ierr)
    if (ierr /= 0) stop "REALLOCATE_IM: problem in attempt to allocate memory"
    if (.not. associated(p)) return
    nold=size(p,1)
    mold=size(p,2)
    real_im_res(1:min(nold,n),1:min(mold,m))=p(1:min(nold,n),1:min(mold,m))
    deallocate(p)
  END FUNCTION reallocate_im

! Initialize or reinitialize the random generator state space to 
! vectors of size LENGTH.  The saved variable SEQ is hashed (via 
! calls to the module routine RAN_HASH) to create unique starting 
! seeds, different for each vector component.

  SUBROUTINE ran_init(length)
    integer(kind=k4b), intent(in) :: length
    integer(kind=k4b) :: new,j,hgt
  
    if (length < lenran) return  !Return if enough space is already allocated.
    hgt=hg
  
! The following lines check that kind value K4B is in fact a 32-bit 
! integer with the usual properties that we expect it to have (under 
! negation and wrap-around addition).  If all of these tests are 
! satisfied, then the routines that use this module are portable, even 
! though they go beyond Fortran 90's integer model.

    if (hg /= 2147483647) stop "RAN_INIT: arith assump 1 fails"
    if (hgng >= 0) stop "RAN_INIT: arith assump 2 fails"
    if (hgt+1 /= hgng) stop "RAN_INIT: arith assump 3 fails"
    if (not(hg) >= 0) stop "RAN_INIT: arith assump 4 fails"
    if (not(hgng) < 0) stop "RAN_INIT: arith assump 5 fails"
    if (hg+hgng >= 0) stop "RAN_INIT: arith assump 6 fails"
    if (not(-1_k4b) < 0) stop "RAN_INIT: arith assump 7 fails"
    if (not(0_k4b) >= 0) stop "RAN_INIT: arith assump 8 fails"
    if (not(1_k4b) >= 0) stop "RAN_INIT: arith assump 9 fails"
  
    if (lenran > 0) then		!Reallocate space, or...
      ranseeds => reallocate_im(ranseeds,length,5)
      ranv => reallocate_iv(ranv,length-1)
      new=lenran+1
    else				!allocate space.
      allocate(ranseeds(length,5))
      allocate(ranv(length-1))
      new=1			 	!Initialize index of first allocation.
      amm=nearest(1.0,-1.0)/hgng	!Use of NEAREST is to ensure that 
    					!returned random #s are < 1.
      if (amm*hgng >= 1.0 .or. amm*hgng <= 0.0) &
       			            stop "RAN_INIT: art assump 10 fails"
    endif

! Set starting values, unique by SEQ and vector component.
    ranseeds(new:,1)=seq
    ranseeds(new:,2:5)=spread(arth(new,1,size(ranseeds(new:,1))),2,4)

! Hash them.
    do j=1,4
      call ran_hash(ranseeds(new:,j),ranseeds(new:,j+1))
    enddo

! Enforce non-negativity.
    where(ranseeds(new:,1:3) <0) ranseeds(new:,1:3)=not(ranseeds(new:,1:3))

! Enforce nonzero.
    where(ranseeds(new:,4:5) == 0) ranseeds(new:,4:5)=1

! Set scalar seeds.
    if (new == 1) then
      iran0=ranseeds(1,1)
      jran0=ranseeds(1,2)
      kran0=ranseeds(1,3)
      mran0=ranseeds(1,4)
      nran0=ranseeds(1,5)
      rans=nran0
    endif

! Point to vector seeds.
    if (length > 1) then
      iran => ranseeds(2:,1)
      jran => ranseeds(2:,2)
      kran => ranseeds(2:,3)
      mran => ranseeds(2:,4)
      nran => ranseeds(2:,5)
      ranv = nran
    endif
    lenran=length
  END SUBROUTINE ran_init

! User interface to release the workspace used by the random 
! number routines

  SUBROUTINE ran_deallocate()
    if (lenran > 0) then
      ranseeds => null()
      ranv => null()
      iran => null()
      jran => null()
      kran => null()
      mran => null()
      nran => null()
      deallocate(ranseeds,ranv)
!      nullify(ranseeds,ranv,iran,jran,kran,mran,nran)
      lenran = 0
    endif
  END SUBROUTINE ran_deallocate

! User interface for seeding the random number routines.  Syntax 
! is exactly like Fortran 90's RANDOM_SEED routine, with one additional 
! argument keywork:  SEQUENCE, set to any integer value, causes an 
! immediate new initialization, seeded by that integer.

  SUBROUTINE ran_seed(sequence,size,put,get)
    integer, optional, intent(in) :: sequence
    integer, optional, intent(out) :: size
    integer, dimension(:), optional, intent(in) :: put
    integer, dimension(:), optional, intent(out) :: get

    if (present(size)) then
      size=5*lenran
    else if (present(put)) then
      if (lenran == 0) return
      ranseeds=reshape(put,shape(ranseeds))

! Enforce non-negativity and nonzero conditions on any user-supplied 
! seeds.
      where (ranseeds(:,1:3) < 0) ranseeds(:,1:3)=not(ranseeds(:,1:3))
      iran0=ranseeds(1,1)
      jran0=ranseeds(1,2)
      kran0=ranseeds(1,3)
      mran0=ranseeds(1,4)
      nran0=ranseeds(1,5)
    else if (present(get)) then
      if (lenran == 0) return
      ranseeds(1,1:5)=(/ iran0,jran0,kran0,mran0,nran0 /)
      get=reshape(ranseeds,shape(get))
    else if (present(sequence)) then
      call ran_deallocate()
      seq=sequence
    endif
  END SUBROUTINE ran_seed

! Hashing of two 32-bit integers, using shifts, xor's, and adds to 
! make the internal nonlinear function.

  SUBROUTINE ran_hash(il,ir)
    integer(kind=k4b), dimension(:), intent(inout) :: il,ir
    integer(kind=k4b), dimension(size(il)) :: is
    integer(kind=k4b) :: j
  
! The various constants are chosen to give good bit mixing and 
! should not be changed.
    do j=1,4
      is=ir
      ir=ieor(ir,ishft(ir,5))+1422217823
      ir=ieor(ir,ishft(ir,-16))+1842055030
      ir=ieor(ir,ishft(ir,9))+80567781
      ir=ieor(il,ir)
      il=is
    enddo
  END SUBROUTINE ran_hash

! RAN2_V from Numerical Recipes in Fortran 90.
!
! Lagged Fibonacci generator combined with a Marsaglia shift 
! sequence and a linear congruential generator.  Returns as HARVEST 
! uniform random numbers between 0.0 and 1.0 (exclusive of the 
! endpoint values).  This generator has the same calling and 
! initialization conventions as Fortran 90's RAND)M_NUIBMER routine.  
! Use RANS_SEED to initialize of reinitialize to a particular sequence.  
! The period of this generator is about 8.5 x 10^37.

  SUBROUTINE ran(harvest)
    real(kind=sp), dimension(:), intent(out) :: harvest
    integer(kind=k4b) :: n

    n=size(harvest)

! Initialization routine in RAN_STATE.
    if (lenran < n+1) call ran_init(n+1)

! Update Fibonacci generator, which has period 4.6 x 10^18.
    ranv(1:n)=iran(1:n)-kran(1:n)
    where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
    iran(1:n)=jran(1:n)
    jran(1:n)=kran(1:n)
    kran(1:n)=ranv(1:n)

! Update Marsaglia shift sequence with period 4.3 x 10^9.
    nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
    nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
    nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))

! Update using shifts instead of multiples.  Wrap-around 
! addition tested at initialization.
    ranv(1:n)=iand(mran(1:n),65535)
    mran(1:n)=ishft(353*ishft(mran(1:n),-16)+ranv(1:n),16) + &
              3533*ranv(1:n)+820265819_k4b
          
! Combine the generators.
    ranv(1:n)=ieor(nran(1:n),kran(1:n))+mran(1:n)

! Make the result positive definite (note that AMM is negative).
    harvest=amm*merge(ranv(1:n),not(ranv(1:n)),ranv(1:n) < 0)
  END SUBROUTINE ran

END MODULE random

