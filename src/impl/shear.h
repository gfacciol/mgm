#ifndef SHEAR_H
#define SHEAR_H

// Horizontal shear and translation of an image
void image_shear(float *in, float *out, int width, int height, float shear,
		 float translation);

#endif				/* !SHEAR_H */
