/* Copyright 2016, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr> */
//
// Created by Gabriele Facciolo on 16/08/16.
//
#include "img_interp.h"
#include "img_tools.h"
#include <cmath>

void leftright_test(struct Img &dx, struct Img &Rdx)
{
    int nc = dx.ncol;
    int nr = dx.nrow;
    int Rnc = Rdx.ncol;
    int Rnr = Rdx.nrow;

    for(int y=0;y<nr;y++)
        for(int x=0;x<nc;x++) {
            int i=x+y*nc;
            int Lx,Rx;
            Lx = round(x+ dx[i]);
            if( (Lx)<Rnc && (Lx)>=0 ){
                int Lidx = Lx + y*Rnc;
                float Rx = Lx + Rdx[Lidx];
                if ( fabs(Rx-x) > 1) {
                    dx[i]  = NAN;
                }
            }else {
                dx[i]  = NAN;
            }
        }

}

void leftright_test_bleyer(struct Img &dx, struct Img &Rdx)
// warps the pixels of the right image to the left, if no pixel in the
// left image receives a contribution then it is marked as occluded
{
    int nc = dx.ncol;
    int nr = dx.nrow;
    int Rnc = Rdx.ncol;
    int Rnr = Rdx.nrow;

    struct Img occL(nc,nr);
    for(int i=0;i<nr*nc;i++) occL[i]=0;

    for(int y=0;y<Rnr;y++)
        for(int x=0;x<Rnc;x++) {
            int i=x+y*Rnc;
            int Lx = round(x+ Rdx[i]);
            if( (Lx)<nc && (Lx)>=0 ){
                occL[Lx + y*nc] = 255;
            }
        }

    for(int i=0;i<nr*nc;i++)
        if(occL[i]==0) dx[i] = NAN;

}

/* min-filter the input with correlation corr, and generates a new
 * disparity map with correlation newcorr.
 * The size of the window is given by w.
 * nan or inf pixels are not considered by the filter.*/
void minfilter(struct Img &disp, struct Img &corr, int w)
{
    int nc = disp.nx;
    int nr = disp.ny;
    /* in case the window size is not odd, we must adapt the right limit.
     * This is because the correlation treats the even windows as if their center O
     * is displaced towards right:  xxxOxx */
    int wl=w/2;
    int wr=w/2;
    if(w%2==0) wr--;

    struct Img outdisp(disp), outcorr(corr);

    for(int y=wl;y<nr-wr;y++)
        for(int x=wl;x<nc-wr;x++) {
            float currcorr=INFINITY;
            for(int i=-wl;i<=wr;i++)
                for(int j=-wl;j<=wr;j++) {
                    if ( currcorr > corr(x+i,y+j) && std::isfinite(disp(x+i,y+j))) {
                        currcorr     = corr(x+i, y+j);
                        outcorr(x,y) = currcorr;
                        outdisp(x,y) = disp(x+i, y+j);
                    }
                }
        }
    disp = outdisp;
    corr = outcorr;
}

void mindiff(struct Img &disp, struct Img &corr, int w, float tau=1.0)
{
    int nc = disp.nx;
    int nr = disp.ny;
    /* in case the window size is not odd, we must adapt the right limit.
     * This is because the correlation treats the even windows as if their center O
     * is displaced towards right:  xxxOxx */
    int wl=w/2;
    int wr=w/2;
    if(w%2==0) wr--;

    struct Img outdisp(disp), outcorr(corr);

    for(int y=wl;y<nr-wr;y++)
        for(int x=wl;x<nc-wr;x++) {
            float mincorr=INFINITY;
            float mindisp=0;
            for(int i=-wl;i<=wr;i++)
                for(int j=-wl;j<=wr;j++) {
                    if(mincorr > corr(x+i, y+j) && std::isfinite(disp(x+i, y+j))) {
                        mincorr = corr(x+i, y+j);
                        mindisp = disp(x+i, y+j);
                    }
                }
            if(std::abs(outdisp(x,y) - mindisp) > tau) {
                outcorr(x,y) = INFINITY;
                outdisp(x,y) = NAN;
            }
        }
    disp = outdisp;
    corr = outcorr;
}


#define SKIP_MAIN
#include "remove_small_cc.c"



// use the current disparity (disp) to update the disparity ranges:  dminI, dmaxI  
// the invalid pixels of disp are assigned the previous disparity range dminP dmaxP
std::pair<float, float> update_dmin_dmax(struct Img disp, struct Img *dminI, struct Img *dmaxI, struct Img &dminP, struct Img &dmaxP, int slack=3, int radius=2) {
    struct Img tdminI(*dminI);
    struct Img tdmaxI(*dmaxI);
    int nx = disp.nx;
    int ny = disp.ny;

    // global (finite) min and max (legacy value for vminP, and vmaxP)
    std::pair<float,float>gminmax = image_minmax(disp);
    float gmin = gminmax.first; float gmax = gminmax.second;

    if (slack<0) slack = -slack;
    int r=radius;

    for (int j=0;j<ny;j++)
        for (int i=0;i<nx;i++)
        {
            float dmin = INFINITY; float dmax = -INFINITY;
            for (int dj=-r;dj<=r;dj++)
                for (int di=-r;di<=r;di++)
                {
                    float v = valneumann(disp, i+di, j+dj);
                    float vminP = valneumann(dminP, i+di, j+dj);
                    float vmaxP = valneumann(dmaxP, i+di, j+dj);
                    if (std::isfinite(v)) {
                        dmin = fmin( dmin, v - slack );
                        dmax = fmax( dmax, v + slack );
                    } else {
                        dmin = fmin( dmin, vminP );
                        dmax = fmax( dmax, vmaxP );
                    }
                }
            if (std::isfinite(dmin)) {
                tdminI[i+j*nx] = dmin; 
                tdmaxI[i+j*nx] = dmax;
            }

        }

    *dminI = tdminI;
    *dmaxI = tdmaxI;
    return std::pair<float, float> (gmin, gmax);
}


// generate the backprojected image u according to the flow
// the image v is used when the flow is invalid at a given point
// when relative_offsets==false the values of flow are absolute image proitions
struct Img backproject_image(struct Img &u, struct Img &v, struct Img &flow, bool relative_offset) {
    struct Img syn = Img(u.nx, u.ny, u.nch);
    for (int x = 0; x < u.nx; x++)
        for (int y = 0; y < u.ny; y++) {
            Point p(x, y);
            if(relative_offset)
                p = p + Point(flow(x,y,0), flow.nch >1 ? flow(x,y,1) : 0 );
            else
                p =     Point(flow(x,y,0), flow.nch >1 ? flow(x,y,1) : y );
            for (int c = 0; c < u.nch; c++)
                if (check_inside_image(p, v))
                    syn.data[x + y * u.nx + c * u.npix] = v.getpixel(p.x, p.y, c);
                else
                    syn.data[x + y * u.nx + c * u.npix] = u.getpixel(x,y,c);
        }
    return syn;
}


