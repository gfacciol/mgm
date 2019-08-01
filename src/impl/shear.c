// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright 2012, Carlo De Franchis <carlo.de-franchis@ens-cachan.fr>
// All rights reserved.

#include <math.h>
#include <fftw3.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

/*
   This function applies shear and translation to the mono-channel image stored at
   address 'in', and write the result at address 'out'. The width and height of
   input image are provided in parameters 'width' and 'height'.

   The transformed image is computed line by line. Each line of the input image
   is translated of a different amout. The basic flow of the algorithm is :
   For each line, compute its DCT, change the phase of each coefficient (this
   results in a tricky formula involving cosinus, sinus and the DST because we use
   DCT and not DFT), then compute the inverse DFT and DST to get the result.
   */

void image_shear(float *in, float *out, int width, int height, float shear, float translation)
{
   // Length of the 1d signal
   int n = width;

   // fftw objects
   double *line_in, *dct, *dst, *line_out_sym, *line_out_antisym;
   fftw_plan p_fwd_dct, p_back_dct, p_back_dst;

   // Memory allocation
   line_in = fftw_malloc(sizeof(double) * n);
   dct = fftw_malloc(sizeof(double) * n);
   dst = fftw_malloc(sizeof(double) * n);
   line_out_sym = fftw_malloc(sizeof(double) * n);
   line_out_antisym = fftw_malloc(sizeof(double) * n);
   if (NULL == line_in || NULL == dct || NULL == dst ||
       NULL == line_out_sym || NULL == line_out_antisym)
      fprintf(stderr, "fftw_malloc() failed in image_shear()\n");

   // Plans computation
   p_fwd_dct = fftw_plan_r2r_1d(n, line_in, dct, FFTW_REDFT10, FFTW_ESTIMATE);
   p_back_dct = fftw_plan_r2r_1d(n, dct, line_out_sym, FFTW_REDFT01, FFTW_ESTIMATE);
   p_back_dst = fftw_plan_r2r_1d(n, dst, line_out_antisym, FFTW_RODFT01, FFTW_ESTIMATE);

   // Compute the output image row by row :
   // each row of the image is translated of a different amount (because of shear)
   for (int row = 0; row < height; row++)
   {
      // Copy the current row
      for (int i = 0; i < n; i++)
         line_in[i] = in[row * n + i];

      // Compute its DCT
      fftw_execute(p_fwd_dct);

      // Normalize the DCT
      for (int i = 0; i < n; i++)
         dct[i] /= n;

      // Compute the modified DCT and DST
      double t = row * shear + translation;
      double a = (M_PI / n) * t;

      // the dct and dst vector are shifted by 1 element,
      // keep the loop computing the same cos and sin
      // to allow cpu optimization
      dct[0] *= cos(0 * a);
      for (int k = 1; k <= (n - 1); k++) {
         dst[k - 1] = dct[k];
         dst[k - 1] *= sin(k * a);
         dct[k] *= cos(k * a);
      }
      // Process separately the last term, because of shifting
      dst[n - 1] = 0;

      // Compute the iDCT and iDST
      fftw_execute(p_back_dct);
      fftw_execute(p_back_dst);

      // Divide by 2 and copy the translated row in the output image
      for (int i = 0; i < n; i++)
         out[row * n + i] = 0.5 * (line_out_sym[i] + line_out_antisym[i]);
   }

   // Free memory
   fftw_destroy_plan(p_fwd_dct);
   fftw_destroy_plan(p_back_dct);
   fftw_destroy_plan(p_back_dst);
   fftw_free(line_in);
   fftw_free(dct);
   fftw_free(dst);
   fftw_free(line_out_sym);
   fftw_free(line_out_antisym);
}
