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
    
    LMAT_MX(0, m, sca_, double)
    LMAT_MX(1, n, sca_, double)
    
    const index_t m_ = static_cast<index_t>(m);
    const index_t n_ = static_cast<index_t>(n);
    
    default_rand_stream& rs = *get_default_rstream();
    
    LMAT_MX_OUT(0, r, marray::numeric_matrix<int32_t>(m_, n_), mat_, int32_t)    
    
    if (nrhs == 3)
    {
        LMAT_MX(2, b, sca_, double)
        int32_t b_ = static_cast<int32_t>(b);
        check_arg(b_ > 0, "The value of b must be positive.");
        
        uniform_int_distr<int32_t> distr(1, b_);
        r = rand_mat(distr, rs, m_, n_);
    }
    else
    {
        LMAT_MX(2, a, sca_, double)
        LMAT_MX(3, b, sca_, double)
        int32_t a_ = static_cast<int32_t>(a);
        int32_t b_ = static_cast<int32_t>(b);
        check_arg(a_ < b_, "The value of a must be less than b.");
        
        uniform_int_distr<int32_t> distr(a_, b_);
        r = rand_mat(distr, rs, m_, n_);
    }
}

