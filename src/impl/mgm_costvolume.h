/* Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>,
 *                     Carlo de Franchis <carlo.de-franchis@ens-cachan.fr>,
 *                     Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>*/
#ifndef COSTVOLUME_H_
#define COSTVOLUME_H_

#include "point.h"
#include <cstring>
#include "img_interp.h"

#include "img_tools.h"

//#include "census_tools.cc"

#define __max(a,b)  (((a) > (b)) ? (a) : (b))
#define __min(a,b)  (((a) < (b)) ? (a) : (b))

// the type of a cost function
typedef float (*cost_t)(Point,Point,const struct Img&,const struct Img&);


inline float computeC_AD_sub( Point p, Point q, const struct Img &u, const struct Img &v) {
	if( !check_inside_image(p,u) ) return INFINITY;
	if( !check_inside_image(q,v) ) return INFINITY;
	float tmp = 0;
	for(int t=0;t<u.nch;t++) {
      float x = u.getpixel(p.x,p.y,t) - v.getpixel(q.x,q.y,t);
      x = __max(x,-x);
      tmp += x ;
   }
	return tmp;
}

inline float computeC_AD( Point p, Point q, const struct Img &u, const struct Img &v) {
	if( !check_inside_image(p,u) ) return INFINITY;
	if( !check_inside_image(q,v) ) return INFINITY;
	float tmp = 0;
	for(int t=0;t<u.nch;t++) {
      float x = val(u,p,t) - val(v,q,t);
      x = __max(x,-x);
      tmp += x ;
   }
	return tmp;
}
inline float computeC_SD( Point p, Point q, const struct Img &u, const struct Img &v) {
   if( !check_inside_image(p,u) ) return INFINITY;
   if( !check_inside_image(q,v) ) return INFINITY;
   float tmp = 0;
   for(int t=0;t<u.nch;t++) {
      float x = val(u,p,t) - val(v,q,t);
      x = __max(x,-x);
      tmp += x * x;
   }
   return tmp;
}


// SMART_PARAMETER(LAMBDA,0.9);

// inline float computeC_COMBI(Point p, Point q, const struct Img &u, const struct Img &v)
// {
// 	float a = computeC_AD(p, q, u, v);
// 	float b = computeC_SC(p, q, u, v);
// 	float l = LAMBDA();
// 	float r = (1 - l) * a + l * b;
// 	return r;
// }


// this macro defines the function CENSUS_NCC_WIN() that reads
// the environment variable with the same name
SMART_PARAMETER(CENSUS_NCC_WIN,3)

float compute_census_distance_array(uint8_t *a, uint8_t *b, int n);

// fast census (input images must pre-processed by census transform)
inline float computeC_census_on_preprocessed_images( Point p, Point q, const struct Img &u, const struct Img &v)
{
	float r = 0;
	for (int t = 0; t < u.nch; t++)
	{
		float vp = valzero(u, p, t);
		float vq = valzero(v, q, t);

	   r += compute_census_distance_array((uint8_t*)(&vp), (uint8_t*)(&vq), sizeof(float));
	}
   // magic factor each channel contrinutes to r  4 bytes (32bits) 
   // to guarantee that the total cost is below 256 we take r*8 / nch
//	return r * 1.0 / u.nch; // magic factor
//	Normalize the costs to 5x5 windows
   const float ratio = 5*5 / (CENSUS_NCC_WIN()*CENSUS_NCC_WIN());
	return r * 1.0 * ratio / u.nch; // magic factor
}


// birchfield and tomasi absolute differences
inline float BTAD( Point p, Point q, int channel, const struct Img &u, const struct Img &v) {
      if( !check_inside_image(p,u) ) return INFINITY;
      if( !check_inside_image(q,v) ) return INFINITY;
      #define min3(a,b,c) (((a)<(b))? (((a)<(c))?(a):(c)) : (((c)<(b))?(c):(b)) )
      #define max3(a,b,c) (((a)>(b))? (((a)>(c))?(a):(c)) : (((c)>(b))?(c):(b)) )

      int t=channel;
      float IL = val(u,p,t);
      float ILp = IL, ILm = IL;
      if (p.x<u.nx-1) ILp = (IL + val(u,p+Point(1,0),t))/2.0;
      if (p.x >= 1  ) ILm = (IL + val(u,p+Point(-1,0),t))/2.0;

      float IR = val(v,q,t);
      float IRp = IR, IRm = IR;
      if (q.x<v.nx-1) IRp = (IR + val(v,q+Point(1,0),t))/2.0;
      if (q.x >= 1  ) IRm = (IR + val(v,q+Point(-1,0),t))/2.0;

      float IminR = min3(IRm,IRp,IR);
      float ImaxR = max3(IRm,IRp,IR);

      float IminL = min3(ILm,ILp,IL);
      float ImaxL = max3(ILm,ILp,IL);

      float dLR =  max3( 0, IL - ImaxR, IminR - IL);
      float dRL =  max3( 0, IR - ImaxL, IminL - IR);

      float BT = __min(dLR, dRL);
      return fabs(BT);
}


// birchfield and tomasi absolute differences
inline float computeC_BTAD( Point p, Point q, const struct Img &u, const struct Img &v) {
	if( !check_inside_image(p,u) ) return INFINITY;
	if( !check_inside_image(q,v) ) return INFINITY;
	float val = 0;
	for(int t=0;t<u.nch;t++)  {
      val += BTAD(p,q,t,u,v);
   }
	return val;
}
// birchfield and tomasi squared differences
inline float computeC_BTSD( Point p, Point q, const struct Img &u, const struct Img &v) {
   if( !check_inside_image(p,u) ) return INFINITY;
   if( !check_inside_image(q,v) ) return INFINITY;
   float val = 0;
   for(int t=0;t<u.nch;t++)  {
      float x = BTAD(p,q,t,u,v);
      val += x*x;
   }
   return val;
}

// Clipped NCC
// NOTE: window size = 3x3
inline float computeC_clippedNCC( Point p, Point q, const struct Img &u, const struct Img &v)
{
   int hwindow = CENSUS_NCC_WIN()/2;
	float r = 0;
   float NCC = 0;
	for (int t = 0; t < u.nch; t++)
	{
      float mu1  = 0; float mu2  = 0;
      float s1   = 0; float s2   = 0;
      float prod = 0;
      int n = 0;
		for (int i = -hwindow; i <= hwindow; i++)
		for (int j = -hwindow; j <= hwindow; j++)
		{
			float v1 = valnan(u, p + Point(i, j), t);
			float v2 = valnan(v, q + Point(i, j), t);
         if (std::isnan(v1) || std::isnan(v2)) return INFINITY;
         mu1+=v1;    mu2+=v2;
         s1 +=v1*v1; s2 +=v2*v2; prod+=v1*v2;
         n++;
      }
      mu1/=n; mu2/=n;
      s1 /=n;  s2/=n; prod/=n;

      NCC += (prod - mu1*mu2) / sqrt( __max(0.0000001,(s1 - mu1*mu1)*(s2 - mu2*mu2)) );
	}
   float clippedNCC = u.nch - __max(0,__min(NCC,u.nch));
	return clippedNCC*64;
}



//// global table of all the cost functions
static struct distance_functions{
   cost_t f;
   const char *name;
} global_table_of_distance_functions[] = {
         #define REGISTER_FUNCTIONN(x,xn) {x, xn}
         REGISTER_FUNCTIONN(computeC_AD,"ad"),
         REGISTER_FUNCTIONN(computeC_SD,"sd"),
         REGISTER_FUNCTIONN(computeC_census_on_preprocessed_images,"census"),
         REGISTER_FUNCTIONN(computeC_clippedNCC,"ncc"),
         REGISTER_FUNCTIONN(computeC_BTAD,"btad"),
         REGISTER_FUNCTIONN(computeC_BTSD,"btsd"),
         REGISTER_FUNCTIONN(computeC_AD_sub,"ad_sub"),
         #undef REGISTER_FUNCTIONN
         {NULL, ""},
};
inline int get_distance_index(const char *name) {
   int r=0; // default cost function is computeC_AD (first in table distance_functions)
   for(int i=0; global_table_of_distance_functions[i].f; i++)
      if (strcmp (name,global_table_of_distance_functions[i].name)==0) 
         r=i;
   return r;
}


//// global table of the prefilter names
static const char* global_table_of_prefilters[] = {
                                             "none",
                                             "census", 
                                             "sobelx",
                                             "gblur",
                                              NULL,
                                           };
inline int get_prefilter_index(const char *name) {
   int r=0;
   for(int i=0; global_table_of_prefilters[i]; i++)
      if (strcmp (name,global_table_of_prefilters[i])==0) 
         r=i;
   return r;
}



#include "dvec.cc"


struct costvolume_t {
   std::vector< Dvec > vectors;
   inline Dvec  operator[](int i) const  { return this->vectors[i];}
   inline Dvec& operator[](int i)        { return this->vectors[i];}
};


// this buffer has the spape of a cost volume
// but can only hold the last maxsize Dvec elements
// inserted with the method set
struct costvolume_buffer_t {
   private:
   std::vector<int> index; // buffer lookup index
   std::vector< Dvec > q;  // circular queue of elements
   int lastq;              // head of the queue, we don't care about the tail

   public:
   int maxsize;            // size of the queue

   costvolume_buffer_t(struct costvolume_t &cv, int m) {
      index   = std::vector<int> (cv.vectors.size()) ;
      q       = std::vector< Dvec > (m);
      lastq   = m-1;
      maxsize = m;
   }
   inline Dvec& get (int i) {
      return q[index[i]];
   }
   inline void set(int i, Dvec& x) {
      lastq = (lastq+1)%maxsize;
      index[i] = lastq;
      q[lastq] = x;
   }
};


struct costvolume_t allocate_costvolume (struct Img min, struct Img max);


struct costvolume_t allocate_and_fill_sgm_costvolume (struct Img &in_u, // source (reference) image                      
                                                      struct Img &in_v, // destination (match) image                  
                                                      struct Img &dminI,// per pixel max&min disparity
                                                      struct Img &dmaxI,
                                                      char* prefilter,        // none, sobel, census(WxW)
                                                      char* distance,         // census, l1, l2, ncc(WxW), btl1, btl2
                                                      float truncDist,        // truncated differences
                                                      float ZOOMFACTOR=1.0);   // subpixel factor (dmin & dmax are stretched)


// computes the second local minimum in the cost volume wrt the fist minimum srtored in dl
// the value of the first minimum is stored in the first channel of cl this function writes its second channel
// THIS IS A LEGACY FUNCTION TO BE REMOVED SOON
void second_local_minimum_from_costvolume(struct costvolume_t &S, struct Img &dl, struct Img *cl);


// computes the PKR (peak ratio) confidence over the cost volume S, for the disparity given in dl
// PKR is the ratio of the second and the first local minimum of the cost
struct Img compute_PKR_confidence(struct costvolume_t &S, struct Img &dl);

// dumps the costvolume in a file with format:
// nx(int x1), ny(int x1), ndisp(int x1), minimum_disp(int x1), costs(float x nx*ny*ndisp)
void dump_costvolume(struct costvolume_t CC, int nx, int ny, int dmin, int dmax, char* cvfilename);
// reads the costvolume from a file with format: 
// nx(int x1), ny(int x1), ndisp(int x1), minimum_disp(int x1), costs(float x nx*ny*ndisp)
void read_costvolume(char* cvfilename, struct costvolume_t &CC, int nx, int ny,  struct Img &dmin, struct Img &dmax);
// computes the right costvolume rearranging the costs contained in the left one
struct costvolume_t right_costvolume_from_left(struct costvolume_t &CCL, int nxL, int nyL, int nxR, int nyR, struct Img &dminR, struct Img &dmaxR);

#endif //COSTVOLUME_H_
