/**********************************************************
 *
 *  logsumexp.cpp
 *
 *  Entropy of discrete distributions
 *
 **********************************************************/

#include <light_mat/matlab/matlab_port.h>
#include <dolphin/common/common_calc.h>

using namespace lmat;
using namespace lmat::matlab;
using namespace dolphin;

LMAT_FP_MEX(logsumexp, 0)
{
    LMAT_MX_NARGINCHK(1, 1)
    LMAT_MX_NARGOUTCHK(0, 1)
    
    LMAT_MX(0, a, mat_, T)
    
    exp_terms<T> et;
    
    if (is_vector(a))
    {       
        et.set_logvalues(a);
        plhs[0] = marray::from_scalar(et.logsum());
    }
    else
    {
        const index_t n = a.ncolumns();
        LMAT_MX_OUT(0, r, marray::numeric_matrix<T>(1, n), row_, T);
        
        for (index_t j = 0; j < n; ++j)
        {
            et.set_logvalues(a.column(j));
            r[j] = et.logsum();
        }
    }
}
