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
    
    LMAT_MX_SCA(0, m, index_t)
    LMAT_MX_SCA(1, n, index_t)
    
    default_rand_stream& rs = *get_default_rstream();
    
    LMAT_MX_OUT(0, r, marray::numeric_matrix<int32_t>(m, n), mat_, int32_t)    
    
    if (nrhs == 3)
    {
        LMAT_MX_SCA(2, b, int32_t)
        check_arg(b > 0, "The value of b must be positive.");
        
        uniform_int_distr<int32_t> distr(1, b);
        r = rand_mat(distr, rs, m, n);
    }
    else
    {
        LMAT_MX_SCA(2, a, int32_t)
        LMAT_MX_SCA(3, b, int32_t)
        check_arg(a < b, "The value of a must be less than b.");
        
        uniform_int_distr<int32_t> distr(a, b);
        r = rand_mat(distr, rs, m, n);
    }
}

