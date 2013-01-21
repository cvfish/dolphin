/**********************************************************
 *
 *  dd_entropy.cpp
 *
 *  Entropy of discrete distributions
 *
 **********************************************************/

#include <light_mat/matlab/matlab_port.h>
#include <dolphin/common/common_calc.h>

using namespace lmat;
using namespace lmat::matlab;
using namespace dolphin;


LMAT_FP_MEX(dd_entropy, 0)
{
    LMAT_MX_NARGINCHK(1, 2)
    LMAT_MX_NARGOUTCHK(0, 1)
    
    LMAT_MX(0, a, mat_, T)
    const index_t m = a.nrows();
    const index_t n = a.ncolumns();
    
    if (nrhs == 1)
    {
        if (m == 1 || n == 1) // as a vector
        {
            T v = entropy(a);
            plhs[0] = marray::from_scalar(v);
        }
        else // colwise
        {
            LMAT_MX_OUT(0, r, marray::numeric_matrix<T>(1, n), row_, T)
            colwise_entropy(a, r);
        }        
    }
    else
    {
        LMAT_MX(1, dim, sca_, double)
        int d = static_cast<int>(dim);
        
        check_arg(d == 1 || d == 2, "dim must be either 1 or 2.");
        
        if (d == 1)  // colwise
        {   
            if (n == 1) // as a vector
            {
                T v = entropy(a);
                plhs[0] = marray::from_scalar(v);
            }
            else // colwise
            {
                LMAT_MX_OUT(0, r, marray::numeric_matrix<T>(1, n), row_, T)
                colwise_entropy(a, r);
            }
        }        
        else
        {                        
            if (m == 1) // as a vector
            {
                T v = entropy(a);
                plhs[0] = marray::from_scalar(v);
            }
            else // rowwise
            {
                LMAT_MX_OUT(0, r, marray::numeric_matrix<T>(m, 1), col_, T)
                rowwise_entropy(a, r);
            }
        }
    }
    
}

