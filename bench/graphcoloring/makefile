CC=mpicc
CXX=mpicxx
LINK=mpicxx
CFLAGS=-I.

%.o: %.c
				$(CC) -c -o $@ $< $(CFLAGS)

%.o: %.cpp
				$(CXX) -c -o $@ $< $(CFLAGS)

all: graphcoloring samplesortdisk samplesort seqsort

seqsort: seqsort.o 
			  $(LINK) -o $@ $< 

samplesort: samplesort.o 
			  $(LINK) -o $@ $< 
samplesortdisk: samplesortdisk.o 
			  $(LINK) -o $@ $< 

graphcoloring: graphcoloring.o 
			  $(LINK) -o $@ $< 

clean:
	rm -f *.o *~ seqsort samplesort samplesortdisk graphcoloring