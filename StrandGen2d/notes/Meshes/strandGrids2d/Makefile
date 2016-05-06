#F90      	= ifort
#FFLAGS   	= -O4 -convert big_endian -r8
F90		= gfortran
#F90FLAGS	= -O4 -fdefault-real-8 -fconvert=big-endian -fPIC
F90FLAGS   	= -O4 -fdefault-real-8 -fconvert=little-endian -fPIC


all: airfoil curve line plate plateNASA bump airfoilNASAtest

airfoil: avdefs
	$(F90) $(F90FLAGS) airfoil.f90 avmeshRoutines.f90 sfuns.F -o airfoil

curve: avdefs
	$(F90) $(F90FLAGS) curve.f90 avmeshRoutines.f90 sfuns.F -o curve

line: avdefs
	$(F90) $(F90FLAGS) line.f90 avmeshRoutines.f90 sfuns.F -o line

plate: avdefs
	$(F90) $(F90FLAGS) plate.f90 avmeshRoutines.f90 sfuns.F -o plate

plateNASA: avdefs
	$(F90) $(F90FLAGS) plateNASA.f90 avmeshRoutines.f90 -o plateNASA

bump: avdefs
	$(F90) $(F90FLAGS) bump.f90 avmeshRoutines.f90 -o bump
	
airfoilNASAtest: avdefs
	$(F90) $(F90FLAGS) airfoilNASAtest.f90 avmeshRoutines.f90 -o airfoilNASAtest

avdefs:
	$(F90) $(F90FLAGS) -c avmeshRoutines.f90

clean:
	rm -f *~ airfoil curve line plate plateNASA bump airfoilNASA *.mod
