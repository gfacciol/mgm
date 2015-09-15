/* Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>,
 *                     Carlo de Franchis <carlo.de-franchis@ens-cachan.fr>,
 *                     Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>*/
#ifndef DVEC_H_
#define DVEC_H_
#include <stdio.h>
#include <math.h>

struct Dvec
{
   float *data;
   int min,max;
   float minval; // minimum value cache

   inline int init(int min, int max, float *dataptr)
   {
      assert(min<max);
   	this->min = min;
   	this->max = max;
      this->data   = dataptr; 
      this->minval = INFINITY;  // by default cache is invalid
      return 0;
   }
//   inline int init(int min, int max)
//   {
//      assert(min<max);
//   	this->min = min;
//   	this->max = max;
//      this->data = std::vector<float >(max-min+1,0);
//      this->minval=INFINITY;  // by default cache is invalid
//      return 0;
//   }
//
//
//   Dvec(int min, int max)
//   {
//      init(min,max);
//   }
//
//   Dvec(int Numel) 
//   {
//      init(0,Numel);
//   }
//
//   Dvec() 
//   {
//   }

   inline float get_minvalue() {
      if (minval == INFINITY) 
         for(int o=min;o<=max;o++) 
            if (this->operator[](o) < minval) 
               minval=this->operator[](o);
      return minval;

   }


   inline void set(int i, float value) {
      if (i>=this->min && i<=this->max) 
      { 
         int idx = i-this->min;
         #pragma omp critical
            data[idx]=value;
         minval=INFINITY;  // invalidate  minval cache
      }
   }

   inline void increment(int i, float value) {
      if (i>=this->min && i<=this->max) 
      { 
         int idx = i-this->min;
         #pragma omp critical
            data[idx]+=value;
         minval=INFINITY;  // invalidate minval cache
      }
   }

   inline void set_nolock(int i, float value) {
      if (i>=this->min && i<=this->max) 
      { 
         int idx = i-this->min;
         data[idx]=value;
         minval=INFINITY;  // invalidate  minval cache
      }
   }

   inline void increment_nolock(int i, float value) {
      if (i>=this->min && i<=this->max) 
      { 
         int idx = i-this->min;
         data[idx]+=value;
         minval=INFINITY;  // invalidate minval cache
      }
   }

   inline float operator[](int i)       { 
      if (i>=this->min && i<=this->max) 
         return data[i-this->min]; 
      else return INFINITY;
   }

};

//int main(){
//   struct Dvec a(10);
//   struct Dvec b(10,20);
//   for(int i=0;i<10;i++){
//      a[0]= i;
//      b[10]= i;
//   }
//   printf("%f\n", a[0]);
//   printf("%f\n", a[13]);
//   printf("%f\n", b[15]);
//}
#endif //DVEC_H_
