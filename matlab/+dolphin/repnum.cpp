/**********************************************************
 *
 *  repnum.cpp
 *
 *  Generate a vector of repeated numbers
 *
 **********************************************************/

#include <light_mat/matlab/matlab_port.h>
#include <light_mat/mateval/mat_reduce.h>

using namespace lmat;
using namespace lmat::matlab;

LMAT_NUM_MEX(repnum, 0)
{
    LMAT_MX_NARGINCHK(1, 2)
    LMAT_MX_NARGOUTCHK(0, 1)
    
    if (nrhs == 1)
    {
        LMAT_MX(0, cnts, mat_, T)
        check_arg(is_vector(cnts), "cnts must be a vector");
        
        const index_t n = cnts.nelems();
        const index_t len = static_cast<index_t>(sum(cnts));
        
        if (is_column(cnts))
        {
            plhs[0] = marray::numeric_matrix<T>(len, 1);
        }
        else
        {
            plhs[0] = marray::numeric_matrix<T>(1, len);
        }
        
        auto r = view_as_col<T>(marray(plhs[0]));
        index_t i = 0;
        for (index_t k = 0; k < n; ++k)
        {
            const index_t cn = static_cast<index_t>(cnts[k]);
            T v = (T)(k + 1);
            
            for (index_t j = 0; j < cn; ++j)
            {
                r[i++] = v;
            }
        }        
    }
    else if (nrhs == 2)
    {
        LMAT_MX(0, x, mat_, T)
        LMAT_MX(1, cnts, mat_, double)
        check_arg(is_vector(x), "x must be a vector");
        check_arg(is_vector(cnts), "cnts must be a vector");
        check_arg(x.nelems() == cnts.nelems(), 
                "x and cnts must have same sizes.");
        
        const index_t n = cnts.nelems();
        const index_t len = static_cast<index_t>(sum(cnts));
        
        if (is_column(cnts))
        {
            plhs[0] = marray::numeric_matrix<T>(len, 1);
        }
        else
        {
            plhs[0] = marray::numeric_matrix<T>(1, len);
        }
        
        auto r = view_as_col<T>(marray(plhs[0]));
        index_t i = 0;
        for (index_t k = 0; k < n; ++k)
        {
            const index_t cn = static_cast<index_t>(cnts[k]);
            T v = x[k];
            
            for (index_t j = 0; j < cn; ++j)
            {
                r[i++] = v;
            }
        } 
    }
}


