#F90      	= ifort
#FFLAGS   	= -O4 -r8 -convert little_endian
F90		= gfortran
F90FLAGS   	= -O4 -fdefault-real-8 -fconvert=little-endian -fPIC


all: lline

lline:
	$(F90) $(F90FLAGS) line.f90 -c -o line.o
	g++ line.o -o line -L/usr/local/gfortran/lib -lgfortran -L$(LD_LIBRARY_PATH)
	$(F90) $(F90FLAGS) sine.f90 -c -o sine.o
	g++ sine.o -o sine -L/usr/local/gfortran/lib -lgfortran -L$(LD_LIBRARY_PATH)
	$(F90) $(FFLAGS) plate.f90 -c -o plate.o
	$(F90) $(FFLAGS) sfuns.F -c -o sfuns.o
	g++ plate.o sfuns.o -o plate -L/usr/local/gfortran/lib -lgfortran -L$(LD_LIBRARY_PATH)
	$(F90) $(F90FLAGS) circle.f90 -c -o circle.o
	g++ circle.o -o circle -L/usr/local/gfortran/lib -lgfortran -L$(LD_LIBRARY_PATH)
	$(F90) $(FFLAGS) bump.f90 sfuns.F -o bump -L$(LD_LIBRARY_PATH)

clean:
	rm -f *~ line plate circle *.o
