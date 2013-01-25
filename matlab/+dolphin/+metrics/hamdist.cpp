/**********************************************************
 *
 *  hamdist.cpp
 *
 *  Hamming distances between vectors
 *
 **********************************************************/

#include "dist_base.h"

LMAT_NUM_AND_BOOL_MEX(hamdist, 0)
{
    LMAT_MX(0, a, mat_, T)
    LMAT_MX(1, b, mat_, T)
    
    distance_port(plhs, dolphin::hamming_distance<T>(), a, b);            
}
