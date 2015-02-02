# Attention
# =========
#
# This makefile exists for the sole purpose of building the software on a particular
# machine, a department's power7. On any other machine, please use cmake, that is
# more convenient.

CC = gcc

# These are the flags to use MPI
MPI_FLAGS = -I/usr/lpp/ppe.poe/include -lmpi -L/usr/lpp/ppe.poe/lib

# The optimization flags to use
OPTIMIZATION_FLAGS = -DNDEBUG -O2

CFLAGS = -Wall -std=c99 $(OPTIMIZATION_FLAGS) $(MPI_FLAGS)

default: gdem

clean:
	rm -r *.o gdem
.PHONY: clean

%.o : %.c
	$(CC) $(CFLAGS) -c -o $@ $<

gdem : parallel.o hll_counter.o parser.o graph.o hyper_anf_mpi.o
	$(CC) $(CFLAGS) -o gdem parallel.o hll_counter.o parser.o graph.o hyper_anf_mpi.o


