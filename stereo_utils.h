//
// Created by Gabriele Facciolo on 16/08/16.
//

#ifndef PROJECT_STEREO_UTILS_CPP_H
#define PROJECT_STEREO_UTILS_CPP_H
//
// Created by Gabriele Facciolo on 16/08/16.
//
#include "img.h"

void leftright_test(struct Img &dx, struct Img &Rdx);

void leftright_test_bleyer(struct Img &dx, struct Img &Rdx);

/* min-filter the input with correlation corr, and generates a new
 * disparity map with correlation newcorr.
 * The size of the window is given by w.
 * nan or inf pixels are not considered by the filter.*/
void minfilter(struct Img &disp, struct Img &corr, int w);

void mindiff(struct Img &disp, struct Img &corr, int w, float tau=1.0);

//#include "remove_small_cc.c"
int remove_small_cc(int w, int h, float *in, float *out, int minarea, float intensity_threshold);

std::pair<float, float> update_dmin_dmax(struct Img outoff, struct Img *dminI, struct Img *dmaxI, struct Img &dminP, struct Img &dmaxP, int slack=3, int radius=2);


struct Img backproject_image(struct Img &u, struct Img &v, struct Img &flow, bool relative_offset=true);
#endif //PROJECT_STEREO_UTILS_CPP_H
