/**********************************************************
 *
 *  rand_gamma.cpp
 *
 *  Generate a matrix of gamma distributed real values
 *
 **********************************************************/

#include "rand_base.h"
#include <light_mat/random/rand_expr.h>

using namespace lmat;
using namespace lmat::matlab;
using namespace lmat::random;

LMAT_SIMPLE_MEX(rand_gamma)
{
    LMAT_MX_NARGINCHK(3, 4)
    LMAT_MX_NARGOUTCHK(0, 1)
    
    LMAT_MX_SCA(0, m, index_t)
    LMAT_MX_SCA(1, n, index_t)
    
    default_rand_stream& rs = *get_default_rstream();
    
    LMAT_MX_OUT(0, r, marray::double_matrix(m, n), mat_, double)    
    
    if (nrhs == 3)
    {
        LMAT_MX_SCA(2, alpha, double)
        
        r = randg(rs, m, n, alpha);
    }
    else 
    {
        LMAT_MX_SCA(2, alpha, double)
        LMAT_MX_SCA(3, beta, double)
        
        r = randg(rs, m, n, alpha, beta);
    }
}

