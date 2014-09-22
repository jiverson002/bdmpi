MPICC = mpicc

LIBS  = -lm

default: ./spmv-mpi;

./spmv-mpi: 
	$(MPICC) -Wall -Wextra -ansi -pedantic -o $@ spmv_util.c spmv_mpi.c $(LIBS)
		
clean:
	rm -f *.o ./spmv-mpi
