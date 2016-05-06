#FORT           = ifort
FORT            = gfortran
#F90FLAGS       = -fpp -O3 -pad -Vaxlib -w -axW -automatic -r8
F90FLAGS        = -O4 -fdefault-real-8 -fconvert=big-endian


OBJS	= fillcsp.o	\
	  fillpsp.o	\
	  findcspdim.o	\
	  findpspdim.o	\
	  initpointingvec.o	\
	  readinput.o	\
	  sfuns.o	\
	  smoothingiter.o	\
	  strand1dDist.o


default: $(OBJS)
	$(FORT) $(F90FLAGS) strand2nc2d.f90 $(OBJS) -o strand2nc2d


clean:
	rm -f *~ *.o strand2nc2d


# suffix rules

.SUFFIXES:      $(SUFFIXES) .F90 .o .f .F .f90 .C .h

.F90.o:
	$(FORT) $(F90FLAGS) $(INCLUDE) -c -o $*.o $*.F90

.f90.o:
	$(FORT) $(F90FLAGS) $(INCLUDE) -c -o $*.o $*.f90

.F.o:
	$(FORT) $(F90FLAGS) $(INCLUDE) -c -o $*.o $*.F

.f.o:
	$(FORT) $(F90FLAGS) $(INCLUDE) -c -o $*.o $*.f

.C.o:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $*.o $*.C
