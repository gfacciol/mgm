CC=gcc -std=c99 
LDFLAGS=
CFLAGS=-O3 -DNDEBUG -ffast-math
#CFLAGS=-g
CXX=g++
CXXFLAGS=$(CFLAGS) -fopenmp

PROGRAMS=mgm_o

all: $(PROGRAMS)


% : %.c 
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

% : %.cc img.cc point.cc 
	$(CXX) $(CXXFLAGS) -DTEST_MAIN $^ -o $@ $(LDFLAGS)

clean:
	rm -f $(PROGRAMS) input.bin output.bin
