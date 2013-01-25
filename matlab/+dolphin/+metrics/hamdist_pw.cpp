/**********************************************************
 *
 *  hamdist_pw.cpp
 *
 *  Hamming distances between vectors
 *
 **********************************************************/

#include "dist_base.h"

LMAT_NUM_AND_BOOL_MEX(hamdist, 0)
{        
    LMAT_MX_NARGINCHK(2, 3)
    LMAT_MX_NARGOUTCHK(0, 1)
    
    LMAT_MX(0, a, mat_, T)
    LMAT_MX(1, b, mat_, T)
    
    if (nrhs == 2)
    {
        pw_distance_port(plhs, dolphin::hamming_distance<T>(), a, b);
    }
    else
    {                
        if (mxIsDouble(prhs[2]))
        {
            LMAT_MX(2, w, col_, double)
            pw_distance_port(plhs,
                    dolphin::weighted_hamming(dolphin::type_<T>(), w),
                    a, b, w.nelems());
        }
        else if (mxIsSingle(prhs[2]))
        {
            LMAT_MX(2, w, col_, float)
            pw_distance_port(plhs,
                    dolphin::weighted_hamming(dolphin::type_<T>(), w),
                    a, b, w.nelems());
        }
        else
        {
            throw lmat::invalid_argument("w must be of floating-point type.");
        }
    }
}
