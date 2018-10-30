/* Copyright 2016, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr> */
#include "mgm_multiscale.h"

#include "mgm_refine.h"
void subpixel_refinement_sgm(struct costvolume_t &S,       // modifies out and outcost
                             std::vector<float > &out,
                             std::vector<float > &outcost,
                             char *refinement); //none, vfit, parabola, cubic, parabolaOCV
#include "mgm_weights.h"
#include "mgm_print_energy.h"
#include "stereo_utils.h"

#include "smartparameter.h"


void zoom_nn(struct Img &in, struct Img *out, int fx, int fy) {
   int nc=in.nx, nr=in.ny, nch=in.nch;
   int onc=out->nx, onr=out->ny;
   // set default
   for(int y=0;y<onc*onr*nch;y++) (*out)[y]=0;

   // copy 
   for(int y=0;y<onr;y++)
   for(int x=0;x<onc;x++)
   for(int c=0;c<nch;c++)
   {
      int xx=x/fx;
      int yy=y/fy;
      (*out)[c*onr*onc + y*onc + x] = in[c*nr*nc + yy*nc + xx];
   }
}

SMART_PARAMETER(MULTISCALE_MINMAX_UPSAMPLE_RADIUS,4);
SMART_PARAMETER(MULTISCALE_MINMAX_UPSAMPLE_SLACK,8); // old value: 3 (it was too small)

void upsample2x_disp(struct Img sdisp, struct Img &refim, struct Img *dmin, struct Img *dmax) {

   for(int i=0; i<sdisp.npix; i++) sdisp[i]*=2.0; // scale disparities

   struct Img outdmin(sdisp.nx, sdisp.ny);
   struct Img outdmax(sdisp.nx, sdisp.ny);
   std::pair<float, float>gminmax = update_dmin_dmax(sdisp, &outdmin, &outdmax, *dmin, *dmax,
         MULTISCALE_MINMAX_UPSAMPLE_SLACK(),
         MULTISCALE_MINMAX_UPSAMPLE_RADIUS()); // only needed to compute dminmax

   zoom_nn(outdmin, dmin, 2, 2);
   zoom_nn(outdmax, dmax, 2, 2);
}



inline double sq(float a,float b) {
   return (a*a+b*b);
}


struct Img downsample2x(struct Img &u, float sigma) {
   // spatial support of the filter and center
   int spatial_support=10;
   int CX=4;

   // allocate output
   struct Img out(ceil(u.nx/2.0),ceil(u.ny/2.0),u.nch);

   // build the filter (gaussian)
   struct Img g(spatial_support,spatial_support);
   for(int j=0;j<g.ny;j++) 
   for(int i=0;i<g.nx;i++) {
      g[i+j*g.nx]=exp( -sq(i-CX-.5,j-CX-.5)/(2.0*sigma*sigma) );
   }

   // apply the filter 
   for(int c=0;c<out.nch;c++) 
   for(int j=0;j<out.ny;j++) 
   for(int i=0;i<out.nx;i++) {

      float acc=0;
      float norm=0;

      for (int y=0;y<g.ny;y++)
      for (int x=0;x<g.nx;x++) {

         if( check_inside_image ( Point((i*2+x-CX), (j*2+y-CX)), u)) {
            float gg = g[x + y*g.nx];
            acc+= u[(i*2+x-CX) + (j*2+y-CX)*u.nx + u.npix*c]*gg;
            norm+=gg;
         }

      }
      out[i+j*out.nx + out.npix*c] = acc/norm;

   }

   return out;
}

// morphologic downsampling 2x of the image u, if is_max is set it applies max, otherwise it's min
inline struct Img downsample2x_disp(struct Img &u, bool is_max) {
   // allocate output
   struct Img o(ceil(u.nx/2.0),ceil(u.ny/2.0),u.nch);

   // apply the filter 
   for(int c=0;c<u.nch;c++) 
   for(int j=0;j<u.ny;j+=2) 
   for(int i=0;i<u.nx;i+=2) {
      float vmin = INFINITY, vmax = -INFINITY;
      for(int k=0;k<2;k++) for(int l=0;l<2;l++) {
         vmin = fmin(vmin, valnan(u, Point(i+k,j+l), c));
         vmax = fmax(vmax, valnan(u, Point(i+k,j+l), c));
      }
      if(is_max)
         o[i/2+j/2*o.nx + o.npix*c] = vmax/2;
      else
         o[i/2+j/2*o.nx + o.npix*c] = vmin/2;
   }
   return o;
}


SMART_PARAMETER(TSGM,4);
SMART_PARAMETER(TSGM_FIX_OVERCOUNT,1);
SMART_PARAMETER(USE_TRUNCATED_LINEAR_POTENTIALS,0);
//SMART_PARAMETER(TSGM_2LMIN,0);

SMART_PARAMETER(REMOVESMALLCC,0.0)
SMART_PARAMETER(MINDIFF,-1)
SMART_PARAMETER(DUMP_COSTVOLUME,0);


////template<pre_function pre, call_function call, post_function post>
////void recursive_multiscale(struct Img &u, struct Img &v, int numscales,
////                          struct Img &dmin, struct Img &dmax, struct Img &dminR, struct Img &dmaxR,
////                          struct Img &dl, struct Img &cl, struct Img &dr, struct Img &cr);



///********************** COSTVOLUME *****************************/
//
//#include "mgm_costvolume.h"
//struct costvolume_t allocate_and_fill_sgm_costvolume (struct Img &in_u, // source (reference) image
//                                                      struct Img &in_v, // destination (match) image
//                                                      struct Img &dminI,// per pixel max&min disparity
//                                                      struct Img &dmaxI,
//                                                      char* prefilter,        // none, sobel, census(WxW)
//                                                      char* distance,         // census, l1, l2, ncc(WxW), btl1, btl2
//                                                      float truncDist,        // truncated differences
//                                                      float ZOOMFACTOR);   // subpixel factor (dmin & dmax are stretched)
//
///********************** MGM *****************************/
//
//#include "mgm_core.cc"
//struct costvolume_t mgm(struct costvolume_t CC, const struct Img &in_w,
//                        const struct Img &dminI, const struct Img &dmaxI,
//                        struct Img *out, struct Img *outcost,
//                        const float P1, const float P2, const int NDIR, const int MGM,
//                        const int USE_FELZENSZWALB_POTENTIALS, // USE SGM(0) or FELZENSZWALB(1) POTENTIALS
//                        int SGM_FIX_OVERCOUNT);                // fix the overcounting in SGM following (Drory etal. 2014)



void mgm_call(struct Img &u, struct Img &v,   // source (reference) image
              struct Img &dmin, struct Img &dmax, struct Img &dminR, struct Img &dmaxR,
              struct Img &dl, struct Img &cl, struct Img &dr, struct Img &cr, struct mgm_param *param)
{
    char* prefilter = param ? param->prefilter : (char*)"none";
    char* refine    = param ? param->refine    : (char*)"none";
    char* distance  = param ? param->distance  : (char*)"ad";
    float truncDist = param ? param->truncDist : INFINITY;
    float P1        = param ? param->P1        : 8;
    float P2        = param ? param->P2        : 32;
    int   NDIR      = param ? param->NDIR      : 4;
    float aP1       = param ? param->aP1       : 1;
    float aP2       = param ? param->aP2       : 1;
    float aThresh   = param ? param->aThresh   : INFINITY;
    float ZOOMFACTOR= param ? param->ZOOMFACTOR: 1.0;

    // special case: receiving a costvolume as input
    if (param->str_dict.count("inputCostVolume") > 0) {
       param->ZOOMFACTOR = 1; // ignore zoom factor
    }

    // compute weights wl and wr that come inside params->img_dict
    struct Img u_w;
    struct Img v_w;
    if (param->img_dict.count("wl")>0 && param->img_dict.count("wr")>0) {
      u_w = compute_mgm_weights_copyvalue(param->img_dict["wl"], aP2, aThresh); // missing aP1 !! TODO
      v_w = compute_mgm_weights_copyvalue(param->img_dict["wr"], aP2, aThresh);
    } else {
      u_w = compute_mgm_weights(u, aP2, aThresh); // missing aP1 !! TODO
      v_w = compute_mgm_weights(v, aP2, aThresh);
    }

    // adapt regularity parameters of USE_TRUNCATED_LINEAR_POTENTIALS==1 when SUBPIX>1
    if (USE_TRUNCATED_LINEAR_POTENTIALS()==1 && ZOOMFACTOR>1) {
       // this conversion for USE_TRUNCATED_LINEAR_POTENTIALS==1 is exact
       P1/=ZOOMFACTOR;
       P2=P2;
    } else {
       // the conversion for USE_TRUNCATED_LINEAR_POTENTIALS==0  is not exact 
       P1/=ZOOMFACTOR;
       P2=P2;
    }

    for(int i = 0; i < TSGM_ITER(); i++)
    {
        //if(scale==0 && i>0) continue;

        //////// START left-right pass
        struct costvolume_t CC, S;
        struct Img zdmin(dmin), zdmax(dmax);
        param->var_dict["right"] = 0;

        // prepare the costvolume 
        if (param->str_dict.count("inputCostVolume") > 0) {
           // load the costvolume from disk and override dmin, dmax
           read_costvolume((char*)param->str_dict["inputCostVolume"].c_str(), CC, u.nx, u.ny, zdmin, zdmax);
           printf("Reading Costvolume: %d %d %f %f\n", u.nx, u.ny, image_minmax(zdmin).first, image_minmax(zdmax).second);
        } else {
           // scale dmin - dmax by the SUBPIX factor and compute costvolume
           for(int i = 0; i < zdmin.npix; i++) zdmin[i] *= ZOOMFACTOR; 
           for(int i = 0; i < zdmax.npix; i++) zdmax[i] *= ZOOMFACTOR;
           CC = allocate_and_fill_sgm_costvolume (u, v, zdmin, zdmax, prefilter, distance, truncDist, ZOOMFACTOR);
        }


        if (DUMP_COSTVOLUME()) {
           int vdmin = image_minmax(zdmin).first;
           int vdmax = image_minmax(zdmax).second;
           printf("%d %d %d %d\n", u.nx, u.ny, vdmin, vdmax);
           dump_costvolume(CC, u.nx, u.ny, vdmin, vdmax,  (char*) "costvolume_left.dat"); 
        }


        //for(int i = 0; i < TSGM_ITER(); i++)
        {
            S = WITH_MGM2() ?
               mgm(CC, u_w, zdmin, zdmax, &dl, &cl, P1, P2,
                    NDIR, TSGM(), param, USE_TRUNCATED_LINEAR_POTENTIALS(), TSGM_FIX_OVERCOUNT()) :
               mgm_naive_parallelism(CC, u_w, zdmin, zdmax, &dl, &cl, P1, P2,
                    NDIR, TSGM(), param, USE_TRUNCATED_LINEAR_POTENTIALS(), TSGM_FIX_OVERCOUNT()) ;
         print_solution_energy(u, dl.data, CC, P1, P2);
        //    std::pair<float,float>gminmax = update_dmin_dmax(dl, &dmin, &dmax, 3);
        //    remove_nonfinite_values_Img(dmin, gminmax.first);
        //    remove_nonfinite_values_Img(dmax, gminmax.second);
        }
        // extract the second local minimum from S
        param->img_dict["confidence_pkrL"] = compute_PKR_confidence(S, dl);
        second_local_minimum_from_costvolume(S, dl, &cl);
        
        // call subpixel refinement  (modifies out and outcost)
        subpixel_refinement_sgm(S, dl.data, cl.data, refine);

        for(int i = 0; i < dl.npix; i++) dl[i] /= ZOOMFACTOR; // scale the solution back by the SUBPIX factor



        //////// START right-left pass
        struct Img zdminR(dminR), zdmaxR(dmaxR); 
        param->var_dict["right"] = 1;

        // prepare the costvolume 
        if (param->str_dict.count("inputCostVolume") > 0) {
           // load the costvolume from disk and override dmin, dmax
           CC = right_costvolume_from_left(CC, u.nx, u.ny, v.nx, v.ny, zdminR, zdmaxR);
           printf("Get Right Costvolume from Left: %d %d %f %f\n", u.nx, u.ny, image_minmax(zdminR).first, image_minmax(zdmaxR).second);
        } else {
           // scale dmin - dmax by the SUBPIX factor and compute costvolume
           for(int i = 0; i < zdminR.npix; i++) zdminR[i] *= ZOOMFACTOR;
           for(int i = 0; i < zdmaxR.npix; i++) zdmaxR[i] *= ZOOMFACTOR;
           CC = allocate_and_fill_sgm_costvolume (v, u, zdminR, zdmaxR, prefilter, distance, truncDist, ZOOMFACTOR);
           // this is actually faster than re-computing the costvolume
           //CC = right_costvolume_from_left(CC, u.nx, u.ny, v.nx, v.ny, zdminR, zdmaxR);
        }


        if (DUMP_COSTVOLUME()) {
           int vdmin = image_minmax(zdminR).first;
           int vdmax = image_minmax(zdmaxR).second;
           printf("%d %d %d %d\n", v.nx, v.ny, vdmin, vdmax);
           dump_costvolume(CC, v.nx, v.ny, vdmin, vdmax,  (char*) "costvolume_right.dat"); 
        }


        //for(int i = 0; i < TSGM_ITER(); i++)
        {
            S = WITH_MGM2() ? 
               mgm(CC, v_w, zdminR, zdmaxR, &dr, &cr, P1, P2,
                    NDIR, TSGM(), param, USE_TRUNCATED_LINEAR_POTENTIALS(), TSGM_FIX_OVERCOUNT()) :
               mgm_naive_parallelism(CC, v_w, zdminR, zdmaxR, &dr, &cr, P1, P2,
                    NDIR, TSGM(), param, USE_TRUNCATED_LINEAR_POTENTIALS(), TSGM_FIX_OVERCOUNT()) ;
         print_solution_energy(v, dr.data, CC, P1, P2);
        //    std::pair<float,float>gminmax = update_dmin_dmax(dr, &dminR, &dmaxR,3);
        //    remove_nonfinite_values_Img(dminR, gminmax.first);
        //    remove_nonfinite_values_Img(dmaxR, gminmax.second);
        }
        // extract the second local minimum from S
        param->img_dict["confidence_pkrR"] = compute_PKR_confidence(S, dr);
        second_local_minimum_from_costvolume(S, dr, &cr);
        
        // call subpixel refinement  (modifies out and outcost)
        subpixel_refinement_sgm(S, dr.data, cr.data, refine);

        for(int i = 0; i < dr.npix; i++) dr[i] /= ZOOMFACTOR;

    }

    //if(scale==0) {
    if(1) {

       if(MEDIAN()) {
          dl = median_filter(dl,MEDIAN());
          dr = median_filter(dr,MEDIAN());
       }

       // TODO: MINDIFF usual values for MINDIFF() ~ 1
       if(MINDIFF()>=0)
          mindiff(dl, cr, CENSUS_NCC_WIN(), MINDIFF());

       // LRRL
       if(TESTLRRL()) {
          Img tmpL(dl);
          Img tmpR(dr);
          leftright_test(dr, tmpL); // R-L
          leftright_test(dl, tmpR); // L-R
       }

       // REMOVE ISOLATED DISPARITY CONNECTED COMPONENTS
       if( REMOVESMALLCC()>0) {
          Img tmpdl(dl), tmpdr(dr);
          remove_small_cc(dl.nx, dl.ny, &tmpdl[0], &dl[0], REMOVESMALLCC(), 5);
          remove_small_cc(dr.nx, dr.ny, &tmpdr[0], &dr[0], REMOVESMALLCC(), 5);
       }
    }
}




void recursive_multiscale(struct Img &u, struct Img &v,
    struct Img &dmin, struct Img &dmax, struct Img &dminR, struct Img &dmaxR,
    struct Img &dl, struct Img &cl, struct Img &dr, struct Img &cr,
    int numscales, int scale, struct mgm_param *param)
{
//    char* prefilter, char* refine, char* distance, float truncDist,
//    const float P1, const float P2, int NDIR,
//    float aP1=1, float aP2=1, float aThresh=INFINITY, int scale=0
//    )

   float maxdisp = -INFINITY, mindisp = INFINITY;
   for (int i=0;i<dmin.npix;i++) {
      maxdisp = fmax(maxdisp, dmax[i]);
      mindisp = fmin(mindisp, dmin[i]);
   }
   float maxdispR = -INFINITY, mindispR = INFINITY;
   for (int i=0;i<dminR.npix;i++) {
      maxdispR = fmax(maxdispR, dmaxR[i]);
      mindispR = fmin(mindispR, dminR[i]);
   }


   if(fmax(u.nx,u.ny) > 100 && fmin(u.nx,u.ny) > 50 && scale < numscales ) {
      struct Img su = downsample2x(u,0.8);
      struct Img sv = downsample2x(v,0.8);
      struct Img sdmin  = downsample2x_disp(dmin, false); // fixme!
      struct Img sdmax  = downsample2x_disp(dmax, true);
      struct Img sdminR = downsample2x_disp(dminR, false);
      struct Img sdmaxR = downsample2x_disp(dmaxR, true);

      struct Img sdl(sdmin);  struct Img scl(sdmin.nx, sdmin.ny,   2);
      struct Img sdr(sdminR); struct Img scr(sdminR.nx, sdminR.ny, 2);

      struct mgm_param sparam(*param); // copy all the params
      // downsample the regularity weight maps
      if (param->img_dict.count("wl")>0 && param->img_dict.count("wr")>0) {
         sparam.img_dict["wl"] = downsample2x(param->img_dict["wl"],0.8);
         sparam.img_dict["wr"] = downsample2x(param->img_dict["wr"],0.8);
      }


      recursive_multiscale(su, sv, sdmin, sdmax, sdminR, sdmaxR,
                            sdl, scl, sdr, scl, numscales, scale+1, &sparam);

      upsample2x_disp(sdl, u, &dmin,  &dmax);
      upsample2x_disp(sdr, v, &dminR, &dmaxR);
   }

//{
//       char name[200]; sprintf(name, "dmax_%02d%02d.tif", scale,0); // DEBUG
//	      iio_write_vector_split(name, dmax); // DEBUG
//}
//{
//       char name[200]; sprintf(name, "dmin_%02d%02d.tif", scale,0); // DEBUG
//	      iio_write_vector_split(name, dmin); // DEBUG
//}

    printf("\n%d/%d %dx%d\t", scale, numscales,u.nx,u.ny);
    printf("maxrange: %f %f\n", image_minmax(dmin).first, image_minmax(dmax).second);
    mgm_call(u, v,
             dmin, dmax, dminR, dmaxR,
             dl, cl, dr, cr, param);

//{
//       char name[200]; sprintf(name, "/tmp/DEBUGd_%02d%02d.tif", scale,0); // DEBUG
//	      iio_write_vector_split(name, dl); // DEBUG
//       sprintf(name, "/tmp/DEBUGw_%02d%02d.tif", scale,0); // DEBUG
//	      iio_write_vector_split(name, param->img_dict["wl"]); // DEBUG
//}

}



