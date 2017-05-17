//
// Created by Gabriele Facciolo on 27/08/16.
//

#ifndef PROJECT_MGM_MULTISCALE_H
#define PROJECT_MGM_MULTISCALE_H

/********************** MGM *****************************/

#include "mgm_core.h"

/********************** MGM *****************************/


//template<pre_function pre, call_function call, post_function post>
//void recursive_multiscale(struct Img &u, struct Img &v, int numscales,
//                          struct Img &dmin, struct Img &dmax, struct Img &dminR, struct Img &dmaxR,
//                          struct Img &dl, struct Img &cl, struct Img &dr, struct Img &cr);

struct mgm_param {
    char* prefilter;
    char* refine;
    char* distance;
    float truncDist;
    float P1, P2;
    int NDIR;
    float aP1, aP2;
    float aThresh;
    float ZOOMFACTOR;
};

void mgm_call(struct Img &u, struct Img &v,
              struct Img &dmin, struct Img &dmax, struct Img &dminR, struct Img &dmaxR,
              struct Img &dl, struct Img &cl, struct Img &dr, struct Img &cr, void *param=NULL);



void recursive_multiscale(struct Img &u, struct Img &v,
                          struct Img &dmin, struct Img &dmax, struct Img &dminR, struct Img &dmaxR,
                          struct Img &dl, struct Img &cl, struct Img &dr, struct Img &cr,
                          int numscales, int scale, void *param=NULL);

#endif //PROJECT_MGM_MULTISCALE_H
