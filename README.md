# MGM: More Global Matching

## Introduction:

This is the MGM code written by Gabriele Facciolo, Carlo de Franchis, and Enric Meinhardt.

**There are two variants, MGM and MGM_multiscale**

If you want to use this software or results obtained with it, the following paper should be cited within your publication:

```
G. Facciolo and C. de Franchis and E. Meinhardt, "MGM: A Significantly More Global Matching for Stereovision", in British Machine Vision Conference, BMVC 2015.

Site : http://dev.ipol.im/~facciolo/mgm/
Email: gabriele.facciolo@cmla.ens-cachan.fr
```



## Overview of the algorithm:

This C++ code can be used for approximately optimizing MRF energies, defined on the 4- or 8-connected image grids, of the form:

![\sum_p C_p(D_p) + \sum_{pq} w_{pq} V (D_p, D_q)](https://latex.codecogs.com/png.latex?%5Cdpi%7B150%7D%20%5Cbg_white%20%5Clarge%20E%28D%29%20%3D%20%5Csum_p%20C_p%28D_p%29%20&plus;%20%5Csum%7Bpq%7D%20w%7Bpq%7D%20V%20%28D_p%2C%20D_q%29)

Where ![C_p()](https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B150%7D%20%5Cbg_white%20%5Clarge%20C_p) denotes the unary term for the p-th node, the variables ![D_p](https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B150%7D%20%5Cbg_white%20%5Clarge%20D_p) take values in the range [0,L-1], ![C_p()](https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B150%7D%20%5Cbg_white%20%5Clarge%20C_p) is represented as a costvolume of size W x H x L.

V() denotes the distance function used for specifying the pairwise potentials, it can take one of these two forms (with params P1,P2):
     1) SMG's potential (Hirschmuller'08)

                    |  0  if |a - b|==0
          Vh(a,b) = | P1  if |a - b|==1
                    | P2  otherwise
  or
     2) any potential described in (Felzenszwalb-Huttenlocher'06),
        this code implements truncated linear and linear (P2=inf)

          Vl(a,b) =  min(P1*|a - b|, P2).

The edge weights are given by w_{pq}, w can actually be used to adapt the parameters P1 and P2 on the pixel basis as:
                ![V (D_p, D_q, P1(w(p)), P2(w(p)) )), P2(w(p)) )](https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B150%7D%20%5Cbg_white%20%5Clarge%20V%20%28D_p%2C%20D_q%2C%20P1%28w%28p%29%29%2C%20P2%28w%28p%29%29%20%29)

But the current implementation just multiplies the potential. The weights are represented as a stack of 8 images. For a pixel p each image of the stack contain the weight to the corresponding neighboring pixel: West, Est, S, N, (NW, NE, SE, SW).

That is the first image just contains the weights for pixels to the left. For 4-connectivity, only the first 4 images are read, while for 8-connectivity the whole stack of 8 images is used.



## Stereo code

The code in this directory uses MGM for stereo.
The option of the mgm program are shown when called without parameters.



### Simplified code and Matlab interface

The `matlab` subdirectory contains a simplified version of MGM, and a matlab wrapper for solving optimization problems as described above.



### Compilation

Tested in Linux, OSX and OpenBSD using gcc and clang

```shell
# In the root directory,

# Make and go into a build folder
$ mkdir -p build
$ cd build

# Create the makefile and build the project
# (Uses OpenMP)
$ cmake ..
$ make
```

You might need to install `fftw3` as it is a dependency

```shell
$ sudo apt install fftw-dev 
```



### Running The Examples

Run these inside the build directory!

#### **MGM**

The following line runs MGM with 8 traversals (-O 8) with 3-neighbor recursion (TSGM=3, mgm is usually 2, but this is NEW!).

     * The cost is measured in absolute differences of the horizontal
       sobel derivatives (sobel_x),
     * for the regularity it uses the FELZENSZWALB potential V
       (USE_TRUNCATED_LINEAR_POTENTIALS=1),
     * then refines disparities with V_fitting and postprocesses with
       MEDIAN filter.
    
    MEDIAN=1 USE_TRUNCATED_LINEAR_POTENTIALS=1  TSGM=3 ./run_mgm -P2 20000 -P1 4 -r -120 -R 30 -p sobel_x -truncDist 63 -s vfit -O 8 ../data/fountain23-im?.png /tmp/{disp,cost}.tif


The following line runs MGM with 8 traversals (-O 8) with 3-neighbor recursion (TSGM=3, mgm is usually 2, but this is NEW!)

      * the cost is CENSUS on 3x3 neighors,
      * for the regularity uses the FELZENSZWALB potential V
        (USE_TRUNCATED_LINEAR_POTENTIALS=1),
      * then refines disparities with V_fitting and postprocesses
        with MEDIAN filter.
    
    MEDIAN=1 CENSUS_NCC_WIN=3 USE_TRUNCATED_LINEAR_POTENTIALS=1  TSGM=3 ./run_mgm -P2 20000 -P1 2 -r -120 -R 30 -t census -s vfit -O 8 ../data/fountain23-im?.png /tmp/{disp,cost}.tif

The following line runs a similar experiment with a satellite image

    OMP_NUM_THREADS=4 MEDIAN=1 CENSUS_NCC_WIN=5 TSGM=3 ./run_mgm -r -22 -R 19 -s vfit -t census -O 8 ../data/rectified_{ref,sec}.tif /tmp/{disp,cost}.tif



**You may similarly invoke `./run_mgm_multi` to run the multiscale variant.** More information about the multiscale variant for MGM can be found in the docs.



### Visualising the Results

As the results of the disparity matching are 32bit float images, you might want to use [pvflip](<https://github.com/gfacciol/pvflip>) to visualise them!



### Parameters

Some of the parameters that are exposed by the library and their explanations can be found [here](https://docs.google.com/spreadsheets/d/1RWZpUHwzIbBTrZQG7sjQVR7uFMCn35e54tAg5_dAWRI/edit?usp=sharing).

Cheers!