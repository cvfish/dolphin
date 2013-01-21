/**********************************************************
 *
 *  rand_pick.cpp
 *
 *  Randomly pick integers without replacement
 *
 **********************************************************/

#include "rand_base.h"
#include <light_mat/random/sample_wor.h>

using namespace lmat;
using namespace lmat::matlab;
using namespace lmat::random;

LMAT_SIMPLE_MEX(rand_pick)
{
    LMAT_MX_NARGINCHK(2, 2)
    LMAT_MX_NARGOUTCHK(0, 1)
    
    LMAT_MX(0, n, sca_, double)
    LMAT_MX(1, k, sca_, double)
    const index_t n_ = static_cast<index_t>(n);
    const index_t k_ = static_cast<index_t>(k);
    
    check_arg(n_ > 0, "n must be a positive integer.");
    check_arg(k >= 0 && k_ <= n_, "k must be in the range [0, n]");
    
    default_rand_stream& rs = *get_default_rstream();
    
    LMAT_MX_OUT(0, r, marray::numeric_matrix<int32_t>(1, k_), row_, int32_t)
    
    if (k == 0) return;
    
    if (n > 1)
    {
        if (k == 1)
        {
            std_uniform_int_distr<int32_t> distr(n);
            r[0] = distr(rs) + 1;
        }
        else if (k == 2)
        {
            std_uniform_int_distr<int32_t> d1(n);
            std_uniform_int_distr<int32_t> d2(n-1);
            
            r[0] = d1(rs) + 1;
            r[1] = d2(rs) + 1;
            if (r[1] >= r[0]) ++r[1];            
        }
        else // k > 2
        {
            int32_t nthres = 16 * k;
            
            if (n > nthres) // use hash-based method
            {
                past_avoid_rand_enumerator<int32_t> e(rs, n, 2 * k);
                for (int32_t i = 0; i < k_; ++i)
                {
                    r[i] = e.next() + 1;
                }                
            } 
            else // use shuffle-based method
            {
                rand_shuffle_enumerator<int32_t> e(rs, n);
                for (int32_t i = 0; i < k_; ++i)
                {
                    r[i] = e.next() + 1;
                }
            }
        }
    }
    else // n == 1
    {
        r[0] = 1;
    }
    
}
