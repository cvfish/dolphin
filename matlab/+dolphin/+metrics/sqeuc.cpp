/**********************************************************
 *
 *  sqeuc.cpp
 *
 *  Squared Euclidean distances between vectors
 *
 **********************************************************/

#include "dist_base.h"

LMAT_FP_MEX(sqeuc, 0)
{
    LMAT_MX(0, a, mat_, T)
    LMAT_MX(1, b, mat_, T)
    
    distance_port(plhs, dolphin::sqeuclidean_distance<T>(), a, b);            
}

