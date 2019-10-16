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
      cv[i].init(floorf(min[i]), ceilf(max[i]));
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


// computes the second local minimum in the cost volume wrt the fist minimum srtored in dl
// the value of the first minimum is stored in the first channel of cl this function writes its second channel
// THIS IS A LEGACY FUNCTION TO BE REMOVED SOON
void second_local_minimum_from_costvolume(struct costvolume_t &S, struct Img &dl, struct Img *cl) {
   if (cl->nch < 2) {
      printf("second_local_minimum: nothing to do here. no second channel\n");
      return;
   }

   for (int i = 0; i < cl->nx * cl->ny; i++) {
      float currdisp  = round(dl[i]);
      float firstmin  = S[i][currdisp];
      float secondmin = INFINITY;
      for(int o=S[i].min;o<=S[i].max;o++) {
         if ( (std::abs(o-currdisp ) > 2)  &&  (secondmin > S[i][o])) {
            secondmin = S[i][o];
         }
      }
      (*cl)[i + cl->nx*cl->ny] = secondmin;
   }

}


// computes the PKR (peak ratio) confidence over the cost volume S, for the disparity given in dl
// PKR is the ratio of the second and the first local minimum of the cost
struct Img compute_PKR_confidence(struct costvolume_t &S, struct Img &dl) {
   struct Img out(dl.nx, dl.ny);
   for (int i = 0; i < dl.nx * dl.ny; i++) {
      float currdisp  = round(dl[i]);
      float firstmin  = S[i][currdisp];
      float secondmin = INFINITY;
      for(int o=S[i].min;o<=S[i].max;o++) {
         if ( (std::abs(o-currdisp ) > 2)  &&  (secondmin > S[i][o])) {
            secondmin = S[i][o];
         }
      }
      out[i] = secondmin/fmax(firstmin , 0.01);
   }
   return out;
}


// dumps the costvolume in a file with format:
// nx(int x1), ny(int x1), ndisp(int x1), minimum_disp(int x1), costs(float x nx*ny*ndisp)
void dump_costvolume(struct costvolume_t CC, int nx, int ny, int dmin, int dmax, char* cvfilename) {
   FILE* fid = fopen(cvfilename, "wb");
   int size[4] = {nx,ny,dmax-dmin+1,dmin};
   fwrite(size, sizeof(int), 4, fid);

   for(int i=0;i<nx*ny;i++) {
      for(int o=dmin;o<=dmax;o++) {
         float value = CC[i][o];
         fwrite(&value, sizeof(float), 1, fid);
      }
   }
   fclose(fid);
}


// reads the costvolume from a file with format: 
// nx(int x1), ny(int x1), ndisp(int x1), minimum_disp(int x1), costs(float x nx*ny*ndisp)
void read_costvolume(char* cvfilename, struct costvolume_t &CC, int nx, int ny, struct Img &zdmin, struct Img &zdmax) {
   FILE* fid = fopen(cvfilename, "rb");
   int meta[4]; // {nx,ny,dmax-dmin+1,dmin};
   size_t r = fread(meta, sizeof(int), 4, fid);

   if (r != 4 || meta[0] != nx || meta[1] != ny)
	   fprintf(stderr, "Bad costvolume header on file \"%s\"\n",cvfilename);

   // overwrite dmin, dmax
   int dmin = meta[3]; 
   int dmax = dmin+meta[2]-1;

   for(int i = 0; i < zdmin.npix; i++) zdmin[i] = dmin;
   for(int i = 0; i < zdmax.npix; i++) zdmax[i] = dmax;

   // create costvolume
   CC.vectors = std::vector< Dvec >(nx*ny);

   // fill costvolume
   for(int i=0;i<nx*ny;i++) {

      CC[i].init(dmin,dmax);
      for(int o=dmin;o<=dmax;o++) {
         float value=0;
         size_t r = fread(&value, sizeof(float), 1, fid);
	 if (r != 1)
		 fprintf(stderr, "short costvolume \"%s\"!\n", cvfilename);
         CC[i].set(o,value);
      }

   }
   fclose(fid);
}


// computes the right costvolume rearranging the costs contained in the left one
struct costvolume_t right_costvolume_from_left(struct costvolume_t &CCL, int nxL, int nyL, int nxR, int nyR, struct Img &dminR, struct Img &dmaxR) {

   for(int i = 0; i < dminR.npix; i++) dminR[i] = INFINITY;
   for(int i = 0; i < dmaxR.npix; i++) dmaxR[i] = -INFINITY; 

   for(int y=0;y<nyL;y++) {
   for(int x=0;x<nxL;x++) {
      // for a pixel in the right image
      int i=x+y*nxL;
      // update dmin and dmax on the left image
      for(int o=CCL[i].min;o<=CCL[i].max;o++) {
         if (x+o >=0 && x+o < nxR) {
             dminR[x+o+y*nxR] = fmin(dminR[x+o+y*nxR], -o);
             dmaxR[x+o+y*nxR] = fmax(dmaxR[x+o+y*nxR], -o);
         }
      }
   }
   }

   // allocate the left costvolume with the right disparity range
   struct costvolume_t CCR = allocate_costvolume (dminR, dmaxR);

   // fill the costs from the right costvolume
   for(int y=0;y<nyR;y++) {
   for(int x=0;x<nxR;x++) {
      int i=x+y*nxR;
      for(int o=CCR[i].min;o<=CCR[i].max;o++) {
         if (x+o >=0 && x+o < nxL) {
            CCR[i].set(o, CCL[x+o +y*nxL][-o]);
         }

      }
   }
   }
   return CCR;

}

