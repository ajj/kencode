	SUBROUTINE initialise_random(iseed)
c random-number generator
		integer, allocatable :: seed(:)
		integer :: n,itime,i,count
		real :: r,random

		call system_clock(count)
		call random_seed(n)
		allocate(seed(n), source=count+37*[(i,i=0,n-1)])
		call random_seed(put=seed)
		deallocate(seed)

	RETURN
	END
