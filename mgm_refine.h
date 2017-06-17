/* Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>,
 *                     Carlo de Franchis <carlo.de-franchis@ens-cachan.fr>,
 *                     Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>*/
#ifndef MGM_REFINE_H_
#define MGM_REFINE_H_
#include <math.h>

#include "refine.h"

//// EXPOSED FUNCTIONS 

typedef float (*refine_t)(const float*, float*, float*);

//// global table of all the refinement functions
struct refinement_functions{
   refine_t f;
   const char *name;
} global_table_of_refinement_functions[] = {
         #define REGISTER_FUNCTIONN(x,xn) {x, xn}
         REGISTER_FUNCTIONN(NULL,                     "none"),
         REGISTER_FUNCTIONN(VfitMinimum,              "vfit"),
         REGISTER_FUNCTIONN(ParabolafitMinimum,       "parabola"),
         REGISTER_FUNCTIONN(CubicfitMinimum,          "cubic"),
         REGISTER_FUNCTIONN(ParabolafitMinimumOpenCV, "parabolaOCV"),
         #undef REGISTER_FUNCTIONN
         {NULL, ""},
};
int get_refinement_index(const char *name) {
   int r=0; // default refinemnt none
   int N = sizeof(global_table_of_refinement_functions)/sizeof(refinement_functions);
   for(int i=0; i<N; i++)
      if (strcmp (name,global_table_of_refinement_functions[i].name)==0) 
         r=i;
   return r;
}




/**
 * refines the disparities(out) and the costs(outcost) by interpolating
 * the location of the minimum in the costvolume S using several methods:
 * none, vfit, parabola, cubic, parabolaOCV
 * */
void subpixel_refinement_sgm(struct costvolume_t &S,       // modifies out and outcost
                             std::vector<float > &out,
                             std::vector<float > &outcost,
                             char *refinement)   //none, vfit, parabola, cubic, parabolaOCV
{
   int N=out.size();

   int idx = get_refinement_index(refinement);
   refine_t refine = global_table_of_refinement_functions[idx].f;
   
   if(refine) {
      if (TSGM_DEBUG()) 
         printf("refinement: using %s\n", global_table_of_refinement_functions[idx].name);
      // subpixel refinement
      #pragma omp parallel for
      for(int i=0;i<N;i++) {
         float minP=out[i];
         float minL=outcost[i];
   
         // can only interpolate if the neighboring disparities are sampled
         int o = minP;
         if( o-1 >= S[i].min && o+2 <= S[i].max) {
            float v[4] = {S[i][o-1], S[i][o], S[i][o+1], S[i][o+2]};
            float deltaX=0;
            refine(v, &minL, &deltaX);
            minP = o + deltaX;

            // this fixes the asymmetry of the cubic interpolation 
            // reverse the signal and compute the minimum
            float vr[4] = {S[i][o+1], S[i][o], S[i][o-1], S[i][o-2]};
            float deltaXr=0;
            float minLr=out[i];
            refine(vr, &minLr, &deltaXr);
            if (minLr < minL) {
               minP = o - deltaXr;
               minL = minLr;
            }

         }
   
         out[i] = minP;
         outcost[i] = minL;
   
      }
   }
}
#endif //MGM_REFINE_H_
