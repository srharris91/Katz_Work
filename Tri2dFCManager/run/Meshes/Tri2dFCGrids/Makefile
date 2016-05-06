#F90	  = ifort
#FFLAGS   = -O4 -convert big_endian -r8
F90      = gfortran
FFLAGS   = -O4 -fdefault-real-8


all: 
	$(F90) $(FFLAGS) channel.f90 -o channel -L$(LD_LIBRARY_PATH)
	$(F90) $(FFLAGS) equiTriPts.f90 channelCub.f90 -o channelCub -L$(LD_LIBRARY_PATH)
	$(F90) $(FFLAGS) vtk2fc.f90 -o vtk2fc -L$(LD_LIBRARY_PATH)
	$(F90) $(FFLAGS) plate.f90 sfuns.F -o plate -L$(LD_LIBRARY_PATH)
	$(F90) $(FFLAGS) circle.f90 sfuns.F -o circle -L$(LD_LIBRARY_PATH)
	$(F90) $(FFLAGS) bump.f90 sfuns.F -o bump -L$(LD_LIBRARY_PATH)
	$(F90) $(FFLAGS) equiTriPts.f90 plateCub.f90 sfuns.F -o plateCub -L$(LD_LIBRARY_PATH)
