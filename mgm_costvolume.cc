/* Copyright 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr> */
#include "point.h"

#include "img_interp.h"

#include "img_tools.h"

#include "mgm_core.h"
#include "mgm_costvolume.h"

#define __max(a,b)  (((a) > (b)) ? (a) : (b))
#define __min(a,b)  (((a) < (b)) ? (a) : (b))

//#include "common/census_tools.cc"
struct Img census_transform(struct Img &bx, int winradius);

extern "C" {
#include "shear.h"
}

/* shift in1 horizontally by q pixels 
 *  *  * (could a noninteger factor) save the translated image in out */
struct Img shift(const struct Img &in1, float q)
{
   int nc = in1.ncol, nr =in1.nrow, nch = in1.nch;
   int N=nc*nr;

   struct Img out(nc,nr,nch);

   float *ti = (float*) malloc(nc*nr*sizeof(float));
   float *to = (float*) malloc(nc*nr*sizeof(float));

   for (int c=0;c<nch;c++){
      for(int i=0;i<N;i++) ti[i] = in1[i+ c*N];
      image_shear(ti, to, nc, nr, 0., -q);
      for(int i=0;i<N;i++) out[i+ c*N] = to[i];
   }
   free(ti);
   free(to);

   return (out);
}

inline int goodmod(int a, int b)
{
       int r = a % b;
       return r < 0 ? r + b : r;
}

std::vector< struct Img > alloc_prefiltered_fourier_subpix_interp(const struct Img &u, int ord)
{
   std::vector< struct Img > us(ord);
   for(int i=0;i<ord;i++){
      us[i] = shift (u, ((float)i)/ ((float)ord));
      //us.push_back(shift (u, ((float)i)/ ((float)ord)));
//       char name[200]; sprintf(name, "/tmp/vv%02d.tif", i); // DEBUG
//	    iio_write_vector_split(name, us[i]); // DEBUG
   }
   return us;
}


struct costvolume_t allocate_costvolume (struct Img min, struct Img max) 
{
   struct costvolume_t cv;
   cv.vectors = std::vector< Dvec >(min.npix);
   for (int i=0;i< min.npix;i++) {
      cv[i].init(min[i], max[i]);
   }

   return cv;
}

struct costvolume_t allocate_and_fill_sgm_costvolume (struct Img &in_u, // source (reference) image                      
                                                      struct Img &in_v, // destination (match) image                  
                                                      struct Img &dminI,// per pixel max&min disparity
                                                      struct Img &dmaxI,
                                                      char* prefilter,        // none, sobel, census(WxW)
                                                      char* distance,         // census, l1, l2, ncc(WxW), btl1, btl2
                                                      float truncDist,        // truncated differences
                                                      float ZOOMFACTOR)   // subpixel factor (dmin & dmax are stretched)
{

   int nx = in_u.nx; 
   int ny = in_u.ny;
   int nch= in_u.nch;

   struct Img u(in_u);
   struct Img v(in_v);
   std::vector< struct Img > vs = alloc_prefiltered_fourier_subpix_interp(v, ZOOMFACTOR);

   // 0. pick the prefilter and cost functions
   int distance_index  = get_distance_index(distance);
   int prefilter_index = get_prefilter_index(prefilter);
   cost_t cost = global_table_of_distance_functions[distance_index].f;

   // 1. parameter consistency check
   if (distance_index == get_distance_index("census") || prefilter_index == get_prefilter_index("census")) {
       if (TSGM_DEBUG()) printf("costvolume: changing both distance and prefilter to CENSUS\n");
       distance_index  = get_distance_index("census");
       prefilter_index = get_prefilter_index("census");
   }
   if (TSGM_DEBUG()) printf("costvolume: selecting distance  %s\n", global_table_of_distance_functions[distance_index].name);
   if (TSGM_DEBUG()) printf("costvolume: selecting prefilter %s\n", global_table_of_prefilters[prefilter_index]);
   if (TSGM_DEBUG()) printf("costvolume: truncate distances at %f\n", truncDist);

   // 2. apply prefilters if needed
   if (prefilter_index == get_prefilter_index("census")) {
      int winradius = CENSUS_NCC_WIN() / 2;
      if (TSGM_DEBUG()) printf("costvolume: applying census with window of size %d\n", winradius*2+1);
      u = census_transform(in_u, winradius);
      //v = census_transform(in_v, winradius);
      for(int i=0;i<ZOOMFACTOR;i++)
         vs[i] = census_transform(vs[i], winradius);

   }
   if (prefilter_index == get_prefilter_index("sobelx")) {
      if (TSGM_DEBUG()) printf("costvolume: applying sobel filter\n" );
      float sobel_x[] = {-1,0,1, -2,0,2, -1,0,1};
      u = apply_filter(in_u, sobel_x, 3, 3, 1);
      //v = apply_filter(in_v, sobel_x, 3, 3, 1);
      for(int i=0;i<ZOOMFACTOR;i++)
         vs[i] = apply_filter(vs[i], sobel_x, 3, 3, 1);
   }
   if (prefilter_index == get_prefilter_index("gblur")) {
      if (TSGM_DEBUG()) printf("costvolume: applying gblur(s=1) filter\n" );
      u = gblur_truncated(in_u, 1.0);
      //v = gblur_truncated(in_v, 1.0);
      for(int i=0;i<ZOOMFACTOR;i++)
         vs[i] = gblur_truncated(vs[i], 1.0);
   }

   // 3. allocate the cost volume 
   struct costvolume_t CC = allocate_costvolume(dminI, dmaxI);

   // 4. apply it 
   #pragma omp parallel for
   for(int jj=0; jj<ny; jj++) for(int ii=0; ii<nx; ii++)
   {
      int pidx  = (ii + jj*nx);
      int allinvalid = 1;

      for(int o=CC[pidx].min;o<=CC[pidx].max;o++) 
      {
         Point p(ii,jj);      // current point on left image
         //Point q = p + Point(o,0); // other point on right image
         //Point q = p + Point(o/ZOOMFACTOR,0); // other point on right image
         Point q = p + Point(floor(o/ZOOMFACTOR),0);
         int zpl = goodmod(o,int(ZOOMFACTOR));
         // 4.1 compute the cost 
         float e = truncDist * u.nch;
         if (check_inside_image(q, vs[0])) {
            //e = cost(p, q, u, v);
            e = cost(p, q, u, vs[zpl]);
         }
         // 4.2 truncate the cost (if needed)
         e = __min(e, truncDist * u.nch);
         // 4.3 store it in the costvolume
         CC[pidx].set_nolock(o, e); // pragma omp critic is inside set
         if(std::isfinite(e)) allinvalid=0;
      }
      // SAFETY MEASURE: If there are no valid hypotheses for this pixel 
      // (ie all hypotheses fall outside the target image or are invalid in some way)
      // then the cost must be set to 0, for all the available hypotheses
      // Leaving inf would be propagated and invalidate the entire solution 
      if (allinvalid) {
         for(int o=CC[pidx].min;o<=CC[pidx].max;o++) 
         {
            CC[pidx].set_nolock(o, 0); // pragma omp critic is inside set
         }
      }
   }
   return CC;
}

