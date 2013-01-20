/**********************************************************
 *
 *  rand_unif.cpp
 *
 *  Generate a matrix of uniformly distributed real values
 *
 **********************************************************/

#include "rand_base.h"
#include <light_mat/random/rand_expr.h>

using namespace lmat;
using namespace lmat::matlab;
using namespace lmat::random;

LMAT_SIMPLE_MEX(rand_unif)
{
    LMAT_MX_NARGINCHK(2, 4)
    LMAT_MX_NARGOUTCHK(0, 1)
    
    LMAT_MX(0, mf, sca_, double)
    LMAT_MX(1, nf, sca_, double)
    
    const index_t m = static_cast<index_t>(mf);
    const index_t n = static_cast<index_t>(nf);
    
    default_rand_stream& rs = *get_default_rstream();
    
    LMAT_MX_OUT(0, r, marray::double_matrix(m, n), mat_, double)    
    
    if (nrhs == 2)
    {
        r = randu(rs, m, n);
    }
    else if (nrhs == 3)
    {
        double a(0);
        LMAT_MX(2, b, sca_, double)
        
        r = randu(rs, m, n, a, b);
    }
    else 
    {
        LMAT_MX(2, a, sca_, double)
        LMAT_MX(3, b, sca_, double)
        
        r = randu(rs, m, n, a, b);
    }
}

