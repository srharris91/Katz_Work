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
FORT = gfortran

#Flags we want to compile with
# Flags for speed (do not turn it in like this, but use it for your time test)
#FFlags = -O3 -march=native -mfpmath=sse
# Flags for errorchecking...this should be use most of the time.
FFlags = -g -fcheck=all -Wall -Wtabs -Wconversion-extra -finit-real=snan -ffpe-trap=invalid,zero,overflow -fbacktrace

#Name of the executable (DO NOT CHANGE THIS)
EXEC = test.exec

#Order of the files you wish to be compiled (replace .f90 with .o)
# Order matters, uper and lowercase matter...avoid spaces in file names
OBJS = system.o main.o


#This is where you define tasks...the first is the default.
exec: $(OBJS)
	@echo ''
	@echo 'Linking [$(EXEC)] from [$(OBJS)]'
	@echo ''	
	@$(FORT) $(FFlags) $(OBJS) -o $(EXEC)


clean:
	@rm *.o $(EXEC) *.mod *.dat
    
# Sufix Rules - This is a bunch of mumbo jumbo you needent understand
.SUFFIXES:
.SUFFIXES: .f90 .o
.f90.o:
	@echo 'Compiling [$@] from [$<] using [$(FORT)]'
	@$(FORT) $(FFlags) -c -o $*.o $*.f90
