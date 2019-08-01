# Multiscale MGM summary

Multiscale may be marginally slower than single scale, but the results are more 
robust and for slanted surfaces it can have a big impact on the computation time.
Census window of 5x5 give smoother disparity maps with less detailed than 3x3.
Using subpixel duplicates the computing time but improves the computed edges.

Recommended setups for S2P from below: #4 or #3 (if results of #4 are too noisy), 
eventually set SUBPIX=1 if time is a concern.

## Tested configurations
    
    # 1. census 5x5 
    # time: ~8s
    REMOVESMALLCC=25 MINDIFF=1 CENSUS_NCC_WIN=5 SUBPIX=1 time ./mgm_multi -S 0 -s vfit -t census data/rectified_???.tif /tmp/w5_S0_sp1_{d,c,b}.tif 

    # 2. census 5x5 + multiscale 
    # time: ~10s
    # multiscale yields less outliers: a more robust result
    REMOVESMALLCC=25 MINDIFF=1 CENSUS_NCC_WIN=5 SUBPIX=1 time ./mgm_multi -S 3 -s vfit -t census data/rectified_???.tif /tmp/w5_S3_sp1_{d,c,b}.tif

    # 3. census 5x5 + multiscale + subpixel
    # time: ~20s
    # subpixel improves the edges 
    REMOVESMALLCC=25 MINDIFF=1 CENSUS_NCC_WIN=5 SUBPIX=2 time ./mgm_multi -S 3 -s vfit -t census data/rectified_???.tif /tmp/w5_S3_sp2_{d,c,b}.tif

    # 4. census 3x3 + multiscale + subpixel 
    # time: ~20s
    # 3x3 window improves the edges, and recovers small features
    REMOVESMALLCC=25 MINDIFF=1 CENSUS_NCC_WIN=3 SUBPIX=2 time ./mgm_multi -S 3 -s vfit -t census data/rectified_???.tif /tmp/w3_S3_sp2_{d,c,b}.tif


## Important parameters

    -S 4             : 5 coarse-to-fine scales 
    -s vfit          : VFIT subpixel refinement algorithm
    CENSUS_NCC_WIN=3 : controls the size of the CENSUS window: 3x3 or 5x5 are recommended  
    REMOVESMALLCC=25 : remove regions with uniform disparity smaller than 25 pixels
    MINDIFF=1        : remove some adherence pixels on a window CENSUS_NCC_WIN
    SUBPIX=2         : performs a 1/2-pixel disparity refinement from the computed the integer disparities


## Multiscale refinement

The disparity is refined with coarse to fine multiscale algorithm. 
The solution computed at the previous resolution is used to 
determine the search range of the next scale. The default 
parameters for computing these bounds are
    MULTISCALE_MINMAX_UPSAMPLE_RADIUS=4
    MULTISCALE_MINMAX_UPSAMPLE_SLACK=8
they control the radius over which the min/max disparity of each
pixel are estimated, and some extra range.

Multiscale may miss a feature if these ranges are not properly adjusted.


## Remarks about census cost

It is important to note that raw census costs are very quantized. 
For instance, a 3x3 window can take only one of 9 values: [0,..8]. 
Although the costs computed by MGM are aggregated somehow they keep
this quantized character. 
Because of this the interpolation of census costs is a bit tricky. 
VFIT interpolation does a nice job when SUBPIX=1, but with SUBPIX>1
it tends to produce flattened surfaces. With SUBPIX=2 this is barely noticeable.
