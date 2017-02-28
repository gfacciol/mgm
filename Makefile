CFLAGS=-Iiio -O3 -DNDEBUG -ffast-math
LDFLAGS=-lpng -ltiff -ljpeg -lm
#CFLAGS=-Iiio -g
CXXFLAGS=$(CFLAGS)

# The following conditional statement appends "-std=gnu99" to CFLAGS when the
# compiler does not define __STDC_VERSION__.  The idea is that many older
# compilers are able to compile standard C when given that option.
# This hack seems to work for all versions of gcc, clang and icc.
CVERSION = $(shell $(CC) -dM -E - < /dev/null | grep __STDC_VERSION__)
ifeq ($(CVERSION),)
CFLAGS := $(CFLAGS) -std=gnu99
endif

# use OpenMP only if not clang
ifeq ($(shell $(CC) -v 2>&1 | grep -c "clang"), 0)
CFLAGS := $(CFLAGS) -fopenmp
endif


PROGRAMS=mgm

all: iio.o $(PROGRAMS)

iio.o: iio/iio.c
	$(CC) $(CFLAGS) -c $^ -o $@
	$(CC) $(CFLAGS) -o iion iio/iio_test_named.c iio.o $(LDFLAGS)

% : %.c iio.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

% : %.cc img.cc point.cc iio.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -f $(PROGRAMS) iio.o iion

test: all
	MEDIAN=1 CENSUS_NCC_WIN=3 USE_TRUNCATED_LINEAR_POTENTIALS=1  TSGM=3 ./mgm -P2 20000 -P1 2 -r -120 -R 30 -t census -s vfit -O 8 data/fountain23-im?.png /tmp/{disp,cost}.tif
	MEDIAN=1 USE_TRUNCATED_LINEAR_POTENTIALS=1  TSGM=3 ./mgm -P2 20000 -P1 4 -r -120 -R 30 -p sobel_x -truncDist 63 -s vfit -O 8 data/fountain23-im?.png /tmp/{disp,cost}.tif
