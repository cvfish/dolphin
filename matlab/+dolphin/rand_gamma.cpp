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
    
    LMAT_MX(0, m, sca_, double)
    LMAT_MX(1, n, sca_, double)
    
    const index_t m_ = static_cast<index_t>(m);
    const index_t n_ = static_cast<index_t>(n);
    
    default_rand_stream& rs = *get_default_rstream();
    
    LMAT_MX_OUT(0, r, marray::double_matrix(m_, n_), mat_, double)    
    
    if (nrhs == 3)
    {
        LMAT_MX(2, alpha, sca_, double)
        
        r = randg(rs, m_, n_, alpha);
    }
    else 
    {
        LMAT_MX(2, alpha, sca_, double)
        LMAT_MX(3, beta, sca_, double)
        
        r = randg(rs, m_, n_, alpha, beta);
    }
}

