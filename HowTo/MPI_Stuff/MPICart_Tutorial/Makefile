comp = mpicxx
cores = 12
run = mpirun
name = MPI_Cart_create

all: compile run

compile:
	$(comp) $(name).cpp -o $(name).exec

run:
	$(run) -n $(cores) $(name).exec


