/**********************************************************
 *
 *  nrmexp.cpp
 *
 *  Normalized exponential terms (softmax)
 *
 **********************************************************/

#include <light_mat/matlab/matlab_port.h>
#include <dolphin/common/common_calc.h>

using namespace lmat;
using namespace lmat::matlab;
using namespace dolphin;

LMAT_FP_MEX(nrmexp, 0)
{
    LMAT_MX_NARGINCHK(1, 1)
    LMAT_MX_NARGOUTCHK(0, 1)
    
    LMAT_MX(0, a, mat_, T)
    LMAT_MX_OUT(0, r, marray_like(a), mat_, T);
    
    exp_terms<T> et;
    
    if (is_vector(a))
    {   
        et.set_logvalues(a);
        et.normalize_to(r);
    }
    else
    {
        const index_t n = a.ncolumns();
        
        for (index_t j = 0; j < n; ++j)
        {
            et.set_logvalues(a.column(j));
            auto rj = r.column(j);
            et.normalize_to(rj);
        }
    }
}
