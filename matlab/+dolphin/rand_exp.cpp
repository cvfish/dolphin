/**********************************************************
 *
 *  rand_exp.cpp
 *
 *  Generate a matrix of exponentially distributed 
 *  real values
 *
 **********************************************************/

#include "rand_base.h"
#include <light_mat/random/rand_expr.h>

using namespace lmat;
using namespace lmat::matlab;
using namespace lmat::random;

LMAT_SIMPLE_MEX(rand_exp)
{
    LMAT_MX_NARGINCHK(2, 3)
    LMAT_MX_NARGOUTCHK(0, 1)
    
    LMAT_MX_SCA(0, m, index_t)
    LMAT_MX_SCA(1, n, index_t)
    
    default_rand_stream& rs = *get_default_rstream();
    
    LMAT_MX_OUT(0, r, marray::double_matrix(m, n), mat_, double)    
    
    if (nrhs == 2)
    {
        r = rande(rs, m, n);
    }
    else
    {
        LMAT_MX_SCA(2, lambda, double)
        
        r = rande(rs, m, n, lambda);
    }
}

