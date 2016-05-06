name = test_vectors_2
name_mod = vector_module_2.F90
comp = gfortran
flags=-Wall -Wextra -Wconversion -fcheck='all' -O3 -o
# -Wall gives Write out of ALL flags (warnings)
#  -Wextra same as above, but for more warnings
#  -Wconversion
#  fcheck='all' checks all fortran checking (includes out of bounds array indices
#  -O3 optimization
#  -o output file
run = ./

all: compile run

compile:
	$(comp) $(name_mod) $(name).F90 $(flags) $(name).exec

run:
	$(run)$(name).exec


