CFLAGS   ?= -O3 -DNDEBUG -ffast-math -march=native
CXXFLAGS := $(CFLAGS)

CPPFLAGS := -Iiio
LDLIBS   := -lfftw3 -lpng -ltiff -ljpeg -lm


BIN       = mgm_multi  mgm
OBJ       = mgm_core.o  mgm_costvolume.o  mgm_multiscale.o  census_tools.o \
            stereo_utils.o  point.o  shear.o  img.o  iio.o

all       : $(BIN)
iio.o     : iio/iio.c       ; $(CC) $(CFLAGS) $(CPPFLAGS) -c $^ -o $@
%         : main_%.o $(OBJ) ; $(CXX) $(LDFLAGS)              $^ -o $@ $(LDLIBS)
clean     :                 ; $(RM) $(BIN) $(OBJ)


test: all
	MEDIAN=1 CENSUS_NCC_WIN=3 USE_TRUNCATED_LINEAR_POTENTIALS=1  TSGM=3 ./mgm -P2 20000 -P1 2 -r -120 -R 30 -t census -s vfit -O 8 data/fountain23-im?.png /tmp/{disp,cost}.tif
	MEDIAN=1 USE_TRUNCATED_LINEAR_POTENTIALS=1  TSGM=3 ./mgm -P2 20000 -P1 4 -r -120 -R 30 -p sobel_x -truncDist 63 -s vfit -O 8 data/fountain23-im?.png /tmp/{disp,cost}.tif


# hack to add -fopenmp only when the CXX compiler supports it
OMPVER := $(shell $(CXX) -dM -E -fopenmp - 2>/dev/null </dev/null |grep _OPENMP)
ifneq ($(OMPVER),)
LDLIBS := $(LDLIBS) -fopenmp
CXXFLAGS := $(CXXFLAGS) -fopenmp
endif

# hack to add -std=gnu99 to very old compilers
CVERSION := $(shell $(CC) -dM -E - < /dev/null | grep __STDC_VERSION__)
ifeq ($(CVERSION),)
CFLAGS := $(CFLAGS) -std=gnu99
endif

