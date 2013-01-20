/**********************************************************
 *
 *  randi_unif.cpp
 *
 *  Generate a matrix of uniformly distributed integers
 *
 **********************************************************/

#include "rand_base.h"
#include <light_mat/random/rand_expr.h>

using namespace lmat;
using namespace lmat::matlab;
using namespace lmat::random;

LMAT_SIMPLE_MEX(randi_unif)
{
    LMAT_MX_NARGINCHK(3, 4)
    LMAT_MX_NARGOUTCHK(0, 1)
    
    LMAT_MX(0, mf, sca_, double)
    LMAT_MX(1, nf, sca_, double)
    
    const index_t m = static_cast<index_t>(mf);
    const index_t n = static_cast<index_t>(nf);
    
    default_rand_stream& rs = *get_default_rstream();
    
    LMAT_MX_OUT(0, r, marray::numeric_matrix<int32_t>(m, n), mat_, int32_t)    
    
    if (nrhs == 3)
    {
        LMAT_MX(2, bf, sca_, double)
        int32_t b = static_cast<int32_t>(bf);
        check_arg(b > 0, "The value of b must be positive.");
        
        std_uniform_int_distr<int32_t> distr(b);
        r = rand_mat(distr, rs, m, n);
    }
    else
    {
        LMAT_MX(2, af, sca_, double)
        LMAT_MX(3, bf, sca_, double)
        int32_t a = static_cast<int32_t>(af);
        int32_t b = static_cast<int32_t>(bf);
        check_arg(a < b, "The value of a must be less than b.");
        
        uniform_int_distr<int32_t> distr(a, b);
        r = rand_mat(distr, rs, m, n);
    }
}

