import ctypes
import numpy as np
from ctypes import c_int, c_float, c_char_p
from numpy.ctypeslib import ndpointer

def mgm_multi_wrapper(im1, im2, dmin, dmax, P1=8, P2=32, prefilter="none",
                      refine="vfit", distance="census", trunc_dist=np.inf,
                      ndir=8, aP1=1, aP2=1, aThresh=5, scales=3,
                      lib_path='libmgm_multi.so'):
    """
    Args:
        im1, im2 (3D numpy array): numpy arrays with the non-interlaced images
        dmin, dmax (ints): min and max disparities
    """
    lib = ctypes.CDLL('libmgm_multi.so')

    nch1, h1, w1 = im1.shape
    nch2, h2, w2 = im2.shape

    lib.mgm_multi.argtypes = (ndpointer(dtype=c_float, shape=(nch1, h1, w1)), c_int, c_int, c_int,
                              ndpointer(dtype=c_float, shape=(nch2, h2, w2)), c_int, c_int, c_int,
                              c_int, c_int, c_float, c_float,
                              c_char_p, c_char_p, c_char_p,
                              c_float, c_int, c_float, c_float, c_float, c_int)
    lib.mgm_multi.restype = ndpointer(dtype=ctypes.c_float, shape=(h1, w1))

    return lib.mgm_multi(im1.astype(np.float32), w1, h1, nch1,
                         im2.astype(np.float32), w2, h2, nch2,
                         dmin, dmax, P1, P2,
                         prefilter.encode(), refine.encode(), distance.encode(),
                         trunc_dist, ndir, aP1, aP2, aThresh, scales)


if __name__ == '__main__':
    import rasterio
    #u = rasterio.open('../../testoutput/output_pair/tiles/row_0000150_height_350/col_0000150_width_350/pair_1/rectified_ref.tif').read()
    #v = rasterio.open('../../testoutput/output_pair/tiles/row_0000150_height_350/col_0000150_width_350/pair_1/rectified_sec.tif').read()
    u = rasterio.open('data/fountain23-imL.png').read()
    v = rasterio.open('data/fountain23-imR.png').read()
    d = mgm_multi_wrapper(u, v, -120, 30)
    import tifffile
    tifffile.imsave('/tmp/d.tif', d)
