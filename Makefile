CC=gcc -std=c99
LDFLAGS=-lpng -ltiff -ljpeg
CFLAGS=-Iiio -O3 -DNDEBUG -ffast-math -fopenmp
#CFLAGS=-Iiio -g
CXX=g++
CXXFLAGS=$(CFLAGS)

PROGRAMS=mgm

all: iio.o $(PROGRAMS)

iio.o: iio/iio.c
	$(CC) $(CFLAGS)  -c $^ -o $@
	$(CC) $(CFLAGS) -o iion iio/iio_test_named.c iio.o $(LDFLAGS)

% : %.c iio.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

% : %.cc img.cc point.cc iio.o
	$(CXX) $(CXXFLAGS) -DTEST_MAIN $^ -o $@ $(LDFLAGS)

clean:
	rm -f $(PROGRAMS) iio.o iion

test: all
	MEDIAN=1 CENSUS_NCC_WIN=3 USE_TRUNCATED_LINEAR_POTENTIALS=1  TSGM=3 ./mgm -P2 20000 -P1 2 -r -120 -R 30 -t census -s vfit -O 8 data/fountain23-im?.png /tmp/{disp,cost}.tif
	MEDIAN=1 USE_TRUNCATED_LINEAR_POTENTIALS=1  TSGM=3 ./mgm -P2 20000 -P1 4 -r -120 -R 30 -p sobel_x -truncDist 63 -s vfit -O 8 data/fountain23-im?.png /tmp/{disp,cost}.tif
