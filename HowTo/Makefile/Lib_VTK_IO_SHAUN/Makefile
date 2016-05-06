####################################################
#                                                  #
#                  Makefile                        #
#                                                  #
####################################################
# A Make file is a file to help compile large projects.  There are several
# ways to approach this, but this is one of the simplest.  There are open 
# source programs that will generate them for you...cmake being the most popular.

# Make files can have variables and follow rules that we define to compile
# things.  You do not need to understand all of the ins and outs of this
# just go with it.  I'll try to help you understand the inputs below.


#Compiler (always gfortran for us)
#FORT = gfortran
FORT = mpif90

#Flags we want to compile with
# Flags for speed (do not turn it in like this, but use it for your time test)
#FFlags = -O3 -march=native -mfpmath=sse
# Flags for errorchecking...this should be use most of the time.
#FFlags = -g -fcheck=all -Wall -Wtabs -Wconversion-extra -finit-real=snan -ffpe-trap=invalid,zero,overflow -fbacktrace
FcompFlags = -cpp -c -J./mod/ -O0 -C -g -fmax-errors=0 -Wall -Wno-array-temporaries -Warray-bounds -Wcharacter-truncation -Wline-truncation -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure -Wunderflow -Wextra -Wuninitialized -fcheck=all -fmodule-private -ffree-line-length-132 -fimplicit-none -ffpe-trap=invalid,overflow -fbacktrace -fdump-core -finit-real=nan  -O3  -DDEBUG -DMPI2 
Fcomp2Flags = -cpp -c -J./mod/ -O0 -C -g -fmax-errors=0 -Wall -Wno-array-temporaries -Warray-bounds -Wcharacter-truncation -Wline-truncation -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure -Wunderflow -Wextra -Wuninitialized -fcheck=all -ffree-line-length-132 -fimplicit-none -ffpe-trap=invalid,overflow -fbacktrace -fdump-core -finit-real=nan  -O3  -DDEBUG -DMPI2 -I/usr/lib/openmpi/include -I/usr/lib/openmpi/lib
FlinkFlags = -O0 -C -g -fmax-errors=0 -Wall -Wno-array-temporaries -Warray-bounds -Wcharacter-truncation -Wline-truncation -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure -Wunderflow -Wextra -Wuninitialized -fcheck=all -ffree-line-length-132 -fimplicit-none -ffpe-trap=invalid,overflow -fbacktrace -fdump-core -finit-real=nan  -O3  -DDEBUG -DMPI2 -pthread -I/usr/lib/openmpi/lib -L/usr//lib -L/usr/lib/openmpi/lib -lmpi_f90 -lmpi_f77 -lmpi -ldl -lhwloc

#Name of the executable (DO NOT CHANGE THIS)
EXEC = test.exec

#Order of the files you wish to be compiled (replace .f90 with .o)
# Order matters, uper and lowercase matter...avoid spaces in file names
OBJS = IR_Precision.o Lib_Pack_Data.o Lib_Base64.o Lib_VTK_IO_Back_End.o Lib_VTK_IO_INI_XML.o Lib_VTK_IO_FLD_XML.o Lib_VTK_IO_GEO_XML.o Lib_VTK_IO_CON_XML.o Lib_VTK_IO_DAT_VAR_XML.o Lib_VTK_IO_END_XML.o Lib_VTK_IO_VTM_XML.o Lib_VTK_IO_PVTK_XML.o Lib_VTK_IO_PVD_XML.o Lib_VTK_IO_INI.o Lib_VTK_IO_GEO.o Lib_VTK_IO_CON.o Lib_VTK_IO_DAT_VAR.o Lib_VTK_IO_END.o Lib_VTK_IO.o  

#This is where you define tasks...the first is the default.
# all: $(OBJS) special exec
exec: $(OBJS) 
	# @echo ''
	# @echo 'Linking [$(EXEC)] from [$(OBJS)]'
	# @echo ''	
	# @$(FORT) $(FlinkFlags) $(OBJS) set3d_mod.o set3d_mpi.o -o $(EXEC)
# 
# special: 
	# @$(FORT) $(Fcomp2Flags) -c -o set3d_mod.o set3d_mod.f90
	# @$(FORT) $(Fcomp2Flags) -c -o set3d_mpi.o set3d_mpi.f90



clean:
	@rm -f *.o $(EXEC)
	@rm -f ./mod/*
    
# Sufix Rules - This is a bunch of mumbo jumbo you needent understand
.SUFFIXES:
.SUFFIXES: .f90 .o
.f90.o:
	@echo 'Compiling [$@] from [$<] using [$(FORT)]'
	@$(FORT) $(FcompFlags) -c -o $*.o $*.f90
