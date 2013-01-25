/**********************************************************
 *
 *  cosdist_pw.cpp
 *
 *  pairwise Squared Euclidean distances between vectors
 *
 **********************************************************/

#include "dist_base.h"

LMAT_FP_MEX(cosdist, 0)
{        
    LMAT_MX_NARGINCHK(2, 2)
    LMAT_MX_NARGOUTCHK(0, 1)
    
    LMAT_MX(0, a, mat_, T)
    LMAT_MX(1, b, mat_, T)
    
    pw_distance_port(plhs, dolphin::cosine_distance<T>(), a, b);            
}

