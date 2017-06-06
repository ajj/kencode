	REAL FUNCTION random(iseed)
c random-number generator
c gfortran intrinsic procedure: no seed, it is initialised in advance with a
c call to initialise_random (in complicated_guide.for). This may not be required
c on a given platform, but does no harm to be included in general.
	CALL random_number(random)
	RETURN
	END
