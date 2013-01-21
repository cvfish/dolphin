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
    
    LMAT_MX(0, m, sca_, double)
    LMAT_MX(1, n, sca_, double)
    
    const index_t m_ = static_cast<index_t>(m);
    const index_t n_ = static_cast<index_t>(n);
    
    default_rand_stream& rs = *get_default_rstream();
    
    LMAT_MX_OUT(0, r, marray::double_matrix(m_, n_), mat_, double)    
    
    if (nrhs == 2)
    {
        r = rande(rs, m_, n_);
    }
    else
    {
        LMAT_MX(2, lambda, sca_, double)
        
        r = rande(rs, m_, n_, lambda);
    }
}

