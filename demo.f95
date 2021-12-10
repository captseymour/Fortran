program demo

! This program demonstrates the use of Fortran's intrinsic
! U(0,1) random number generator.

implicit none
integer, dimension(8) :: old,seed     ! Needs at least size 8.
integer :: i
real, dimension(3) :: harvest

! To have control of the sequence, we need to 'put' the seed.

seed=12345
call random_seed(put=seed)

! See that the seed is what we set it to.

call random_seed(get=old)
print*, "Old starting value: ",old

! Look at some random uniforms.  There will only be 3 of them.

call random_number(harvest)
print*,"Random numbers: ",harvest

do i=1,3     ! On 3 separate occasions...

! See how the seed changed from the last call of random_number

  call random_seed(get=old)
  print*,"Present starting value: ",old

! We will get a different set of random numbers with that new seed

  call random_number(harvest)
  print*,"Random number: ",harvest

! And since that call changed the seed, we will get yet another 
! set of random numbers.

  call random_number(harvest)
  print*,"Random number: ",harvest
end do

end program demo
