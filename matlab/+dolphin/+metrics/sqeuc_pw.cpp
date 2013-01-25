/**********************************************************
 *
 *  sqeuc_pw.cpp
 *
 *  pairwise Squared Euclidean distances between vectors
 *
 **********************************************************/

#include "dist_base.h"

LMAT_FP_MEX(sqeuc, 0)
{        
    LMAT_MX_NARGINCHK(2, 3)
    LMAT_MX_NARGOUTCHK(0, 1)
    
    LMAT_MX(0, a, mat_, T)
    LMAT_MX(1, b, mat_, T)
    
    if (nrhs == 2)
    {
        pw_distance_port(plhs, dolphin::sqeuclidean_distance<T>(), a, b);
    }
    else
    {
        LMAT_MX(2, w, col_, T)
        pw_distance_port(plhs, dolphin::weighted_sqeuclidean(w), a, b, w.nelems());
    }
}
