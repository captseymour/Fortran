program demo

implicit none
integer, dimension(1) :: old,seed
integer :: i
real, dimension(3) :: harvest

seed(1)=12345
call random_seed(put=seed)
call random_seed(get=old)
print*, "Old starting value: ",old

call random_number(harvest)
print*,"Random numbers: ",harvest

do i=1,3
  call random_seed(get=old)
  print*,"Present starting value: ",old
  call random_number(harvest)
  print*,"Random number: ",harvest
  call random_number(harvest)
  print*,"Random number: ",harvest
end do

end program demo
