/**********************************************************
 *
 *  ldist.cpp
 *
 *  L-p distance between vectors
 *  (subsumes, cityblock/euclidean/chebyshev)
 *
 **********************************************************/

#include "dist_base.h"

LMAT_FP_MEX(ldist, 0)
{
    LMAT_MX_NARGINCHK(2, 4)
    LMAT_MX_NARGOUTCHK(0, 1)
    
    LMAT_MX(0, a, mat_, T)
    LMAT_MX(1, b, mat_, T)
    
    double p_(2);
    if (nrhs > 2)
    {
        LMAT_MX_SCA(2, p, double)
        p_ = p;
    }
    
    if (nrhs <= 3) // non-weighted
    {
        if (p_ == 2)
        {
            distance_port(plhs, dolphin::euclidean_distance<T>(), a, b);
        }
        else if (p_ == 1)
        {
            distance_port(plhs, dolphin::cityblock_distance<T>(), a, b);
        }
        else if (p_ > 0 && math::isinf(p_))
        {
            distance_port(plhs, dolphin::chebyshev_distance<T>(), a, b);
        }
        else
        {
            distance_port(plhs, dolphin::minkowski_distance<T>(T(p_)), a, b);
        }
    }
    else // weighted
    {                
        LMAT_MX(3, w, col_, T)
        const index_t d = w.nelems();
                
        if (p_ == 2)
        {
            distance_port(plhs, dolphin::weighted_euclidean(w), a, b, d);
        }
        else if (p_ == 1)
        {
            distance_port(plhs, dolphin::weighted_cityblock(w), a, b, d);
        }
        else if (p_ > 0 && math::isinf(p_))
        {
            distance_port(plhs, dolphin::chebyshev_distance<T>(), a, b, d);
        }
        else
        {
            distance_port(plhs, dolphin::weighted_minkowski(T(p_), w), a, b, d);
        }
    }

                   
}

