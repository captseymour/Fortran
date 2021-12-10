program demo

use random
implicit none
integer, dimension(1) :: old
integer :: i,k,seed
real, dimension(3) :: harvest

seed=12345
call ran_seed(sequence=seed)
call ran_seed(get=old)
print*, "Old starting value: ",old

call ran(harvest)
print*,"Random numbers: ",harvest

do i=1,3
  call ran_seed(get=old)
  print*,"Present starting value: ",old
  call ran(harvest)
  print*,"Random number: ",harvest
  call ran(harvest)
  print*,"Random number: ",harvest
end do

end program demo

