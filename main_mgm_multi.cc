/* Copyright 2016, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr> */
#include "stdio.h"
#include <vector>
#include <iostream>
#include "assert.h"

#include "smartparameter.h"

//// a structure to wrap images
// // not used here but generally useful
// typedef std::vector<float> FloatVector;
#include "img_interp.h"
#include "img_tools.h"


#include "mgm_multiscale.h"
#include "stereo_utils.h"


SMART_PARAMETER(SUBPIX,1.0)


// c: pointer to original argc
// v: pointer to original argv
// o: option name after hyphen
// d: default value (if NULL, the option takes no argument)
static char *pick_option(int *c, char ***v, char *o, char *d)
{
    int argc = *c;
    char **argv = *v;
    int id = d ? 1 : 0;
    for (int i = 0; i < argc - id; i++)
        if (argv[i][0] == '-' && 0 == strcmp(argv[i] + 1, o)) {
            char *r = argv[i + id] + 1 - id;
            *c -= id + 1;
            for (int j = i; j < argc - id; j++)
                (*v)[j] = (*v)[j + id + 1];
            return r;
        }
    return d;
}




int main(int argc, char* argv[])
{
    /* patameter parsing - parameters*/
    if(argc<4)
    {
        fprintf (stderr, "too few parameters\n");
        fprintf (stderr, "  usage: %s  [-r dmin -R dmax] [-m dminImg -M dmaxImg] [-O NDIR: 2, (4), 8, 16] u v out [cost [backflow]]\n",argv[0]);
        fprintf (stderr, "    [-P1 (8) -P2 (32)]: sgm regularization parameters P1 and P2\n");
        fprintf (stderr, "    [-p PREFILT(none)]: prefilter = {none|census|sobelx|gblur} (census is WxW)\n");
        fprintf (stderr, "    [-t      DIST(ad)]: distance = {census|ad|sd|ncc|btad|btsd}  (ncc is WxW, bt is Birchfield&Tomasi)\n");
        fprintf (stderr, "    [-truncDist (inf)]: truncate distances at nch*truncDist  (default INFINITY)\n");
        fprintf (stderr, "    [-s  SUBPIX(none)]: subpixel refinement = {none|vfit|parabola|cubic}\n");
        fprintf (stderr, "    [-aP1         (1)]: multiplier factors of P1 and P2 when\n");
        fprintf (stderr, "    [-aP2         (1)]:    \\sum |I1 - I2|^2 < nch*aThresh^2\n");
        fprintf (stderr, "    [-aThresh     (5)]: Threshold for the multiplier factor (default 5)\n");
        fprintf (stderr, "    [-S scales    (3)]: Number of scales\n");
        fprintf (stderr, "    [-Rd  fname]: right disparity map\n");
        fprintf (stderr, "    [-Rc  fname]: right cost map\n");
		  fprintf (stderr, "    [-wl  fname]: regularization weights for the left disparity map\n");
		  fprintf (stderr, "    [-wr  fname]: regularization weights for the right disparity map\n");
        fprintf (stderr, "    [-confidence_costL      fnameL -confidence_costR      fnameR]: left and right cost confidence maps\n");
        fprintf (stderr, "    [-confidence_pkrL       fnameL -confidence_pkrR       fnameR]: left and right PKR confidence maps\n");
        fprintf (stderr, "    [-confidence_consensusL fnameL -confidence_consensusR fnameR]: left and right consensus confidence maps\n");
        fprintf (stderr, "    [-inputCostVolume filename]: file containing the costvolume of the left image\n");
        fprintf (stderr, "    ENV: CENSUS_NCC_WIN=3   : size of the window for census and NCC\n");
        fprintf (stderr, "    ENV: TESTLRRL=1   : lrrl\n");
		  fprintf (stderr, "    ENV: REMOVESMALLCC=0 : remove connected components of disp. smaller than (recomended 25)\n");
		  fprintf (stderr, "    ENV: MINDIFF=-1   : remove disp. inconsistent with minfilter on a window of size CENSUS_NCC_WIN (recommended 1)\n");
        fprintf (stderr, "    ENV: TSGM=4       : regularity level\n");
        fprintf (stderr, "    ENV: TSGM_ITER=1  : iterations\n");
        fprintf (stderr, "    ENV: TSGM_FIX_OVERCOUNT=1   : fix overcounting of the data term in the energy\n");
        fprintf (stderr, "    ENV: TSGM_DEBUG=0 : prints debug informtion\n");
//		fprintf (stderr, "      ENV: TSGM_2LMIN=0 : use the improved TSGM cost only for TSGM=2. Overrides TSGM value\n");
        fprintf (stderr, "    ENV: SUBPIX=1     : subpixel steps\n");
        fprintf (stderr, "    ENV: USE_TRUNCATED_LINEAR_POTENTIALS=0 : use the Felzenszwalb-Huttenlocher\n");
        fprintf (stderr, "                      : truncated linear potential (when=1). P1 and P2 change meaning\n");
        fprintf (stderr, "                      : The potential they describe becomes:  V(p,q) = min(P2,  P1*|p-q|)\n");
        return 1;
    }


    //read the parameters
    int i = 1;
    char *in_min_disp_file = pick_option(&argc, &argv, (char*) "m", (char*) "");
    char *in_max_disp_file = pick_option(&argc, &argv, (char*) "M", (char*) "");
    int dmin = atoi(pick_option(&argc, &argv, (char*) "r", (char*) "-30"));
    int dmax = atoi(pick_option(&argc, &argv, (char*) "R", (char*) "30"));
    int NDIR  = atoi(pick_option(&argc, &argv, (char*) "O", (char*) "4"));
    float P1  = atof(pick_option(&argc, &argv, (char*) "P1", (char*) "8"));
    float P2  = atof(pick_option(&argc, &argv, (char*) "P2", (char*) "32"));
    float aP1 = atof(pick_option(&argc, &argv, (char*) "aP1", (char*) "1"));
    float aP2 = atof(pick_option(&argc, &argv, (char*) "aP2", (char*) "1"));
    float aThresh = atof(pick_option(&argc, &argv, (char*) "aThresh", (char*) "5"));

    char* distance  = pick_option(&argc, &argv, (char*) "t", (char*) "ad");   //{census|ad|sd|ncc|btad|btsd}
    char* prefilter = pick_option(&argc, &argv, (char*) "p", (char*) "none"); //{none|census|sobelx}
    char* refine    = pick_option(&argc, &argv, (char*) "s", (char*) "none"); //{none|vfit|parabola|cubic}
    float truncDist = atof(pick_option(&argc, &argv, (char*) "truncDist",  (char*) "inf"));
    int   scales    = atoi(pick_option(&argc, &argv, (char*) "S",  (char*) "3"));

    char* wl_name   = pick_option(&argc, &argv, (char*) "wl", (char*) "");   //weights left
    char* wr_name   = pick_option(&argc, &argv, (char*) "wr", (char*) "");   //weights right
    char* inputCV   = pick_option(&argc, &argv, (char*) "inputCostVolume", (char*)"");

    
    // catch all the other optional output filenames
    char * const keyword_parameters[] = {(char*) "confidence_consensusL", (char*) "confidence_consensusR", // path consensus confidence
                                         (char*) "confidence_costL",      (char*) "confidence_costR",  // cost confidence
                                         (char*) "confidence_pkrL",       (char*) "confidence_pkrR",  // cost confidence
                                         (char*) "Rd",                    (char*) "Rc"}; // right disparity and cost maps
    int num_keyword_parameters = 8;
    std::vector< std::pair< std::string, std::string > >  other_keyword_parameters;
    for(int i = 0; i<num_keyword_parameters; i++) {
       char* fname = pick_option(&argc, &argv, keyword_parameters[i], (char*)"");
       other_keyword_parameters.push_back( std::pair< std::string, std::string>( keyword_parameters[i], fname )  );
    }
    

    char* f_u     = (argc>i) ? argv[i] : NULL;      i++;
    char* f_v     = (argc>i) ? argv[i] : NULL;      i++;
    char* f_out   = (argc>i) ? argv[i] : NULL;      i++;
    char* f_cost  = (argc>i) ? argv[i] : NULL;      i++;
    char* f_back  = (argc>i) ? argv[i] : NULL;      i++;

    printf("%d %d\n", dmin, dmax);


    // read input
    struct Img u = iio_read_vector_split(f_u);
    struct Img v = iio_read_vector_split(f_v);

    remove_nonfinite_values_Img(u, 0);
    remove_nonfinite_values_Img(v, 0);

    struct Img dminI(u.nx, u.ny);
    struct Img dmaxI(u.nx, u.ny);
    for(int i=0;i<u.npix;i++) {dminI[i]=dmin; dmaxI[i]=dmax;}

    if(strcmp (in_min_disp_file,"")!=0 ){
        dminI = iio_read_vector_split(in_min_disp_file);
        dmaxI = iio_read_vector_split(in_max_disp_file);
        // sanity check for nans
        remove_nonfinite_values_Img(dminI, dmin);
        remove_nonfinite_values_Img(dmaxI, dmax);

        // more hacks to prevent errors from bad inputs (min>=max)
        for (int i=0;i<u.npix;i++) {
            if (dmaxI[i] < dminI[i] + 1) dmaxI[i] = ceil(dminI[i] + 1);
        }
    }


    P1 = P1*u.nch; //8
    P2 = P2*u.nch; //32

    // output variables
    struct Img outoff  = Img(u.nx, u.ny);
    struct Img outcost = Img(u.nx, u.ny, 2);

    // variables for LR
    struct Img outoffR  = Img(v.nx, v.ny);
    struct Img outcostR = Img(v.nx, v.ny, 2);
    struct Img dminRI(v.nx, v.ny);
    struct Img dmaxRI(v.nx, v.ny);
    for(int i = 0; i < v.npix; i++) {dminRI[i] = -dmax; dmaxRI[i] = -dmin;}

    // handle multiscale
    struct mgm_param param = {prefilter, refine, distance,truncDist,P1,P2,NDIR,aP1,aP2,aThresh,1.0};
    // load weights for the regularization term
    if(strcmp (wl_name,"")!=0 && strcmp (wr_name,"")!=0) {
       param.img_dict["wl"] = iio_read_vector_split(wl_name);
       param.img_dict["wr"] = iio_read_vector_split(wr_name);
    }
    // input Costvolume
    if(strcmp(inputCV, "") != 0) {
       scales = 0; //disable multiscale if costvolume is given as input
       param.str_dict["inputCostVolume"] = inputCV;
    }

    recursive_multiscale(u,v,dminI,dmaxI,dminRI,dmaxRI,outoff, outcost, outoffR, outcostR, scales, 0, &param);

    // handle subpixel refinement 
    if(SUBPIX()>1) {
       // disparity range is estimated from min,max on a 9x9 window and enlarged by +-2 
       update_dmin_dmax(outoff,  &dminI,  &dmaxI , dminI,  dmaxI,  2, 4);
       update_dmin_dmax(outoffR, &dminRI, &dmaxRI, dminRI, dmaxRI, 2, 4);
       struct mgm_param param = {prefilter, refine, distance,truncDist,P1,P2,NDIR,aP1,aP2,aThresh,(float)SUBPIX()};
       recursive_multiscale(u,v,dminI,dmaxI,dminRI,dmaxRI,outoff, outcost, outoffR, outcostR, 0, 0, &param);
    }

    // save the disparity
    iio_write_vector_split(f_out, outoff);
    // generate the backprojected image
    struct Img syn = backproject_image(u, v, outoff);
    if(f_cost) iio_write_vector_split(f_cost, outcost);
    if(f_back) iio_write_vector_split(f_back, syn);


    param.img_dict["Rd"] = outoffR;
    param.img_dict["Rc"] = outcostR;


    for(int i = 0; i < num_keyword_parameters; i++) {
       std::string keyword =  other_keyword_parameters[i].first;
       std::string fname   =  other_keyword_parameters[i].second;
       if ( fname != "" ) 
          if (param.img_dict.count(keyword)>0) 
             iio_write_vector_split( (char*) fname.c_str(), param.img_dict[keyword.c_str()]);
          else 
             printf("Ignoring keywork %s NOT FOUND in param.img_dict\n", keyword.c_str());
    }


    return 0;
}
