	SUBROUTINE initialise_random(iseed)
c random-number generator

c gfortran intrinsic procedure: no seed, initialised in advance with 
c a call to random_seed (in complicated_guide.for)
	CALL random_number(r)
	random=r

	RETURN
	END
