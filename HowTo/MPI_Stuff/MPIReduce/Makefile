comp = mpif90
cores = 8
run = mpirun
name = reduce_avg

all: compile run

compile:
	$(comp) $(name).f90 -o $(name).exec

run:
	$(run) -n $(cores) $(name).exec


