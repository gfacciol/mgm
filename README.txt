============
Introduction
============
This is the MGM code written by Gabriele Facciolo,
Carlo de Franchis, and Enric Meinhardt.

If you want to use this software or results obtained with it,
the following paper should be cited within your publication:

- G. Facciolo and C. de Franchis and E. Meinhardt,
  "MGM: A Significantly More Global Matching for Stereovision",
  in British Machine Vision Conference, BMVC 2015.

Site : http://dev.ipol.im/~facciolo/mgm/
Email: gabriele.facciolo@cmla.ens-cachan.fr


=========================
Overview of the algorithm
=========================
This C++ code can be used for approximately optimizing MRF energies,
defined on the 4- or 8-connected image grids, of the form:

  E(D) = \sum_p C_p(D_p) + \sum_{pq} w_{pq} V (D_p, D_q)

  where C_p() denotes the unary term for the p-th node, the variables
  D_p take values in the range [0,L-1], C_p is represented as a
  costvolume of size W x H x L.
  V() denotes the distance function used for specifying the pairwise
  potentials, it can take one of these two forms (with params P1,P2):
     1) SMG's potential (Hirschmuller'08)

                    |  0  if |a - b|==0
          Vh(a,b) = | P1  if |a - b|==1
                    | P2  otherwise
  or
     2) any potential described in (Felzenszwalb-Huttenlocher'06),
        this code implements truncated linear and linear (P2=inf)

          Vl(a,b) =  min(P1*|a - b|, P2).

  The edge weights are given by w_{pq}, w can actually be used to adapt
  the parameters P1 and P2 on the pixel basis as:
                V (D_p, D_q, P1(w(p)), P2(w(p)) ),
  but the current implementation just multiplies the potential.
  The weights are represented as a stack of 8 images. For a pixel p
  each image of the stack contain the weight to the corresponding
  neighboring pixel: West, Est, S, N, (NW, NE, SE, SW).
  That is the first image just contains the weights for pixels to the left.
  For 4-connectivity, only the first 4 images are read, while for
  8-connectivity the whole stack of 8 images is used.


===========
Stereo code
===========
  The code in this directory uses MGM for stereo.
  The option of the mgm program are shown when called without parameters.


====================================
Simplified code and Matlab interface
====================================
  The matlab subdirectory contains a simplified version of MGM,
  and a matlab wrapper for solving optimization problems as
  described above.


===========
Compilation
===========
  Run make to compile the C++ code (uses OpenMP)
   Tested in Linux, OSX and OpenBSD using gcc and clang


=============
Example calls
=============

The following line runs MGM with 8 traversals (-O 8) with 3-neighbor
     recursion (TSGM=3, mgm is usually 2, but this is NEW!).

     * The cost is measured in absolute differences of the horizontal
       sobel derivatives (sobel_x),
     * for the regularity it uses the FELZENSZWALB potential V
       (USE_TRUNCATED_LINEAR_POTENTIALS=1),
     * then refines disparities with V_fitting and postprocesses with
       MEDIAN filter.

    MEDIAN=1 USE_TRUNCATED_LINEAR_POTENTIALS=1  TSGM=3 ./mgm -P2 20000 -P1 4 -r -120 -R 30 -p sobel_x \
      -truncDist 63 -s vfit -O 8 data/fountain23-im?.png /tmp/{disp,cost}.tif


The following line runs MGM with 8 traversals (-O 8) with 3-neighbor recursion (TSGM=3, mgm is usually 2, but this is NEW!)

      * the cost is CENSUS on 3x3 neighors,
      * for the regularity uses the FELZENSZWALB potential V
        (USE_TRUNCATED_LINEAR_POTENTIALS=1),
      * then refines disparities with V_fitting and postprocesses
        with MEDIAN filter.

    MEDIAN=1 CENSUS_NCC_WIN=3 USE_TRUNCATED_LINEAR_POTENTIALS=1  TSGM=3 ./mgm -P2 20000 -P1 2 -r -120 -R 30 \
      -t census -s vfit -O 8 data/fountain23-im?.png /tmp/{disp,cost}.tif

The following line runs a similar experiment with a satellite image

    OMP_NUM_THREADS=4 MEDIAN=1 CENSUS_NCC_WIN=5 TSGM=3 ./mgm -r -22 -R 19 -s vfit \
      -t census -O 8 data/rectified_{ref,sec}.tif /tmp/{disp,cost}.tif

====================
Full list of options
====================

This section describes the options accepted by the ``mgm`` program. 

-----
Usage
-----

  mgm <options> left_image.tif right_image.tif output_disparity.tif \
      [cost_function.tif [back_flow.tif]]

--------------------
Command-line options
--------------------

-r (default=-30): 
    Minimum horizontal disparity value. (The images are assumed
    to be rectified, which eliminates the vertical disparity.)

-R (default=30): 
    Maximum horizontal disparity value. 

-O (default=4):
    Number of search directions. Options: 2, 4, 8, 16. 

-P1 (default=8)
    SGM regularization parameter P1.

-P2 (default=32): 
    SGM regularization parameter P2.

-p (default=none): 
    Prefilter algorithm. Options: none, census, sobelx, gblur. The
    ``census`` mode uses the window of dimensions ``CENSUS_NCC_WIN``.

-t (default=ad): 
    Distance function. Options: census, ad, sd, ncc, btad, btsd. For
    ``ncc`` the window of dimensions ``CENSUS_NCC_WIN`` is used. The
    ``bt`` option is the Birchfield-Tomasi distance.

-truncDist (default=inf): 
    Truncate distances at nch * truncDist.

-s (default=none):
    Subpixel refinement algorithm. Options: none, vfit, parabola,
    cubic.

-aP1 (default=1): 
    Multiplier factor of P1 when sum |I1 - I2|^2 < nch * aThresh^2.

-aP2 (default=1): 
    Multiplier factor of P2 as above.

-aThresh (default=5):
   Threshold for the multiplier factors.

-m FILE (default=none): 
    A file with minimum input disparity.

-M FILE (default=none):
    A file with maximum input disparity.
 
-l FILE (default=none): 
    Write here the disparity without the left-to-right test.

-----------------------
Environmental variables
-----------------------

These should be set on the command line before ``mgm`` is invoked.

CENSUS_NCC_WIN=3: 
    Size of the window for the census prefilter algorithm and NCC
    (normalized cross-correlation).

TESTLRRL=1: 
    If 1, do left-to-right and right-to-left consistency checks.

MEDIAN=0:
     Radius of the median filter postprocessing.

TSGM=4:
    Regularity level.

TSGM_ITER=1: 
    Number of iterations.

TSGM_FIX_OVERCOUNT=1: 
    If 1, fix overcounting of the data term in the energy.

TSGM_DEBUG=0:
    If 1, print debug information.

TSGM_2LMIN=0:
    Use the improved TSGM cost only for TSGM=2. Overrides the TSGM
    value.

USE_TRUNCATED_LINEAR_POTENTIALS=0: 
    If 1, use the Felzenszwalb-Huttenlocher truncated linear
    potential. Then P1 and P2 change meaning. The potential they
    describe becomes V(p,q) = min(P2, P1*|p-q|).
