=============
Introduction:
=============
This is a demo code of MGM written by Gabriele Facciolo.

If you want to use this software or results obtained with it, 
the following paper should be cited within your publication:

- G. Facciolo and C. de Franchis and E. Meinhardt,
  "MGM: A Significantly More Global Matching for Stereovision", 
  in British Machine Vision Conference, BMVC 2015.

Site : http://dev.ipol.im/~facciolo/mgm/
Email: gabriele.facciolo@cmla.ens-cachan.fr


=========================
Overview of the C++ code:
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


==============================
Using the C++ code (mgm_o.cc):
==============================
The current C++ implementation is meant for research purposes, so its 
interfaces are not too user friendly. The main function call is:

   mgm(C, ew, dminI, dmaxI, output.data, outcost.data, P1, P2, NDIR, MGM, USE_FELZENSZWALB_V);

where the main parameters are: the cost volume C, the edge weights ew, 
the output labels, and the aggregated costs of the output labels.
The following code describes how to setup all the inputs and extract 
the outputs. Similar calls are used in the standalone wrapper mgm_o.cc.


  // variables TO BE SET
  int W=256, H=256, L=64;     // width, height, #labels
  int NDIR = 4;               // directions 4 or 8
  int MGM  = 2;               // 1 (SGM) or 2 (MGM)
  float P1 = 8, P2 = 32;      // MGM params 
  int USE_FELZENSZWALB_V = 0; // 1 Felzenszwalb trunc linear; 0 SGM's potential

  // create and load the costvolume with YOUR COSTS
  struct Img dmin(W, H), dmax(W, H);
  for ( int i = 0; i < W*H; i++ ) { 
    dmin[i]=0; dmax[i]=L-1;                      // <<- LABEL RANGE
  }
  std::vector< Dvec > C = allocate_costvolume (dminI, dmaxI); 
  for ( int i = 0; i < W*H; i++ ) 
    for ( int o = 0; o < L; o++) 
      C[i].set_nolock(o, _YOURCOSTS[i + o*W*H]); // <<- YOUR COSTS

  // create and load YOUR EDGE WEIGHTS
  struct Img ew(W, H, 8);
  for ( int o = 0; o < 8; o++) 
     for ( int i = 0; i < W*H; i++ ) 
        ew[i + o*W*H] = 1.0;                    // <<- ALL EDGE WEIGHTS = 1.0

  // create output structures
  struct Img output(W, H), outcost(W, H); 

  // call MGM
  mgm(C, ew, dminI, dmaxI, output.data, outcost.data, P1, P2, NDIR, MGM, USE_FELZENSZWALB_V);

  // label readout 
  for ( int i = 0; i < W*H; i++ ) 
     printf("%f\n", output[i]); // <<- RECOVER THE RESULTS


=======================================
Matlab wrapper and example application:
=======================================
The Matlab wrapper MGM_wrapper.m calls the C++ implementation (mgm_o.cc):

  labeling = MGM_wrapper(unary, NDIR, P1, P2, MGM, VTYPE, w)

  INPUTS: 
  unary  the cost volume C of size W x H x L
  NDIR   pass directions 4:
                         8: (default) 
  P1,P2  interaction potential params: P1=8, P2=32 (default)
  MGM    # of messages   1: 1D propagation in each pass as SGM
                         2: use messages from 2 neighbors (default)
  VTYPE  V potential     0: Vh as in SGM (Hirschmuller'08) (default)
                         1: Vl truncated lienar (Felzenszwalb-Huttenlocher'06)
  w      edge weights of size W x H x 8 (default: all ones)

The function stereomatch_MGM.m provides an usage example of the MGM wrapper.
Running the script runme.m, calls stereomatch_MGM with different parameters
and shows the corresponding results.


===========
Compilation
===========
1. run make to compile the C++ code (uses OpenMP)
   Only tested in Linux and OSX using gnu-gcc.

2. call runme.m from Matlab
