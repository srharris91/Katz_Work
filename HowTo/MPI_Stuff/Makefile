comp = mpicxx
cores = 4
run = mpirun
name = MPI_Trial_5

all: compile run

compile:
	$(comp) $(name).cpp -o $(name).exec

run:
	$(run) -n $(cores) $(name).exec


