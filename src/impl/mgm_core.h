#ifndef MGM_H_
#define MGM_H_

#include "smartparameter.h"

SMART_PARAMETER(TSGM_DEBUG,0)
SMART_PARAMETER(TSGM_ITER,1)
SMART_PARAMETER(TESTLRRL,1)

SMART_PARAMETER(MEDIAN,0)
SMART_PARAMETER(WITH_MGM2,0);

#include <map>
#include <string>
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
    std::map<std::string, struct Img> img_dict; // input and output dictionary of images
    std::map<std::string, float> var_dict;
    std::map<std::string, std::string> str_dict;
};

#include "mgm_costvolume.h"

// mgm returns the "aggregated" cost volume, out, and outcost without any other refinement
// This is the diagonally parallel implementation of MGM. Compared to the naive
// version, this one uses less memory and scales better with the number of cores.
// However, I've observed that it may also be more cache intensive.
struct costvolume_t mgm(struct costvolume_t CC, const struct Img &in_w,
                        const struct Img &dminI, const struct Img &dmaxI,
                        struct Img *out, struct Img *outcost,
                        const float P1, const float P2, const int NDIR, const int MGM, struct mgm_param *param,
                        const int USE_FELZENSZWALB_POTENTIALS = 0, // USE SGM(0) or FELZENSZWALB(1) POTENTIALS
                        int SGM_FIX_OVERCOUNT = 1);                 // fix the overcounting in SGM following (Drory etal. 2014)


// mgm returns the "aggregated" cost volume, out, and outcost without any other refinement
// This is the naive parallel implementation of MGM, all traversals (up to 8) are computed
// in parallel. Thus lots of memory is required.
struct costvolume_t mgm_naive_parallelism(struct costvolume_t CC, const struct Img &in_w,
                                          const struct Img &dminI, const struct Img &dmaxI,
                                          struct Img *out, struct Img *outcost,
                                          const float P1, const float P2, const int NDIR, const int MGM, struct mgm_param *param,
                                          const int USE_FELZENSZWALB_POTENTIALS = 0, // USE SGM(0) or FELZENSZWALB(1) POTENTIALS
                                          int SGM_FIX_OVERCOUNT = 1);                 // fix the overcounting in SGM following (Drory etal. 2014)
#endif
