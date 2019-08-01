/* You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.*/
/* Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr> */
#ifndef DVEC_H_
#define DVEC_H_
#include <stdio.h>
#include <math.h>


// this object mimics a vector of floats (on access time at least)
// but with arbitrary index values. The actual storage just for 
// indices in the range [min, max].
//    Dvec v  = Dvec(10,20); // internally just a vector with 11 elements
//    float p = v[11];      // p is 0
//    float q = v[1];       // 1 is out of range so q is INFINITY
struct Dvec
{
   std::vector<float > data;
   int min,max;
   float minval; // minimum value cache

   inline int init(int min, int max)
   {
      assert(min<max);
   	this->min = min;
   	this->max = max;
      this->data = std::vector<float >(max-min+1,0); // devault value is 0 (I know a waste of time but vectors are initialized anyway)
      this->minval = INFINITY;  // by default cache is invalid
      return 0;
   }


   inline Dvec(int min, int max)
   {
      init(min,max);
   }

   inline Dvec(int Numel) 
   {
      init(0,Numel);
   }

   inline Dvec() 
   {
   }

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
//         #pragma omp critical
            data[idx]=value;
         minval=INFINITY;  // invalidate  minval cache
      }
   }

   inline void increment(int i, float value) {
      if (i>=this->min && i<=this->max) 
      { 
         int idx = i-this->min;
//         #pragma omp critical
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

   inline float operator[](int i)       { if (i>=this->min && i<=this->max) return this->data[i-this->min]; else return INFINITY;}

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
