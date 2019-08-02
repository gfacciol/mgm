#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "smartparameter.h"

/* Original QNM's version from E. Meinhardt-Llopis
 * Megawave port by G. Facciolo*/

/*/ a smart parameter is just like a regular parameter, but it can be
// re-defined at the shell-environment.  Instead of
//
// #define NUMBER 42
// ...
// printf("%g", NUMBER);
//
// do
// SMART_PARAMETER(NUMBER,42)
// ...
// printf("%g", NUMBER());
//
// Notice that the environment only gets queried once, at the first use.
//
*/

// For defining global double smart parameters
#define SMART_PARAMETER_GLOBAL(n,v) int smapa_known_ ## n = false;\
double smapa_value_ ## n = v;

// For defining global integer smart parameters
#define SMART_PARAMETER_GLOBAL_INT(n,v) int smapa_known_ ## n = false;\
int smapa_value_ ## n = v;

SMART_PARAMETER_GLOBAL(TSGM_DEBUG,0)
SMART_PARAMETER_GLOBAL(TSGM_ITER,1)
SMART_PARAMETER_GLOBAL(TESTLRRL,1)
SMART_PARAMETER_GLOBAL(MEDIAN,0)
SMART_PARAMETER_GLOBAL(WITH_MGM2,0);
SMART_PARAMETER_GLOBAL(MULTISCALE_MINMAX_UPSAMPLE_RADIUS,4);
SMART_PARAMETER_GLOBAL(MULTISCALE_MINMAX_UPSAMPLE_SLACK,8); // old value: 3 (it was too small)
SMART_PARAMETER_GLOBAL(TSGM,4);
SMART_PARAMETER_GLOBAL(TSGM_FIX_OVERCOUNT,1);
SMART_PARAMETER_GLOBAL(USE_TRUNCATED_LINEAR_POTENTIALS,0);
SMART_PARAMETER_GLOBAL(REMOVESMALLCC,0.0)
SMART_PARAMETER_GLOBAL(MINDIFF,-1)
SMART_PARAMETER_GLOBAL(DUMP_COSTVOLUME,0);
SMART_PARAMETER_GLOBAL(CENSUS_NCC_WIN,3)
SMART_PARAMETER_GLOBAL(SUBPIX,1.0)
SMART_PARAMETER_GLOBAL(TSGM_2LMIN,0);
