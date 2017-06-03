complicated_guide : complicated_guide.o cylindrical_guide.o guide_section.o bin.o cross.o refl.o real_def.o read_char_line.o read_real_line.o makefilename.o initialise_random.o random.o
	gfortran -L/usr/X11/lib -lX11 -Wl,-framework -Wl,Foundation -L/usr/lib -lpng -lz -L/usr/local/pgplot -lpgplot -o complicated_guide complicated_guide.o cylindrical_guide.o guide_section.o bin.o cross.o refl.o real_def.o read_char_line.o read_real_line.o makefilename.o initialise_random.o random.o

complicated_guide.o : complicated_guide.for
	gfortran -c complicated_guide.for

cylindrical_guide.o : cylindrical_guide.for
	gfortran -c cylindrical_guide.for

guide_section.o : guide_section.for
	gfortran -c guide_section.for

bin.o : bin.for
	gfortran -c bin.for

cross.o : cross.for
	gfortran -c cross.for

refl.o : refl.for
	gfortran -c refl.for

real_def.o : real_def.for
	gfortran -c real_def.for

read_char_line.o : read_char_line.for
	gfortran -c read_char_line.for

read_real_line.o : read_real_line.for
	gfortran -c read_real_line.for

makefilename.o : makefilename.for
	gfortran -c makefilename.for

initialise_random.o : initialise_random.for
	gfortran -c initialise_random.for

random.o : random.for
	gfortran -c random.for

