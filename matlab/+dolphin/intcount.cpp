/**********************************************************
 *
 *  intcount.cpp
 *
 *  Counting integer indices/subscripts
 *
 **********************************************************/

#include <light_mat/matlab/matlab_port.h>
#include <dolphin/common/dpaccum.h>

using namespace lmat;
using namespace lmat::matlab;


LMAT_NUM_MEX(intcount, 1)
{
    LMAT_MX_NARGINCHK(2, 3)
    LMAT_MX_NARGOUTCHK(0, 1)
    
    if (nrhs == 2)
    {
        LMAT_MX(0, K, sca_, double)
        LMAT_MX(1, I, col_, T)
        
        const index_t K_ = static_cast<index_t>(K);
        check_arg(K > 0, "K must be a positive integer.");
        
        LMAT_MX_OUT(0, r, marray::numeric_matrix<uint32_t>(K_, 1), 
                col_, uint32_t);
        
        dolphin::add_counts(I - T(1), r);
    }
    else // nrhs == 3
    {
        LMAT_MX(0, siz, row_, double)
        LMAT_MX(1, I, col_, T)
        LMAT_MX(2, J, col_, T)
                
        check_arg(siz.nelems() == 2, "siz must be a pair of numbers.");
        const index_t m = static_cast<index_t>(siz[0]);
        const index_t n = static_cast<index_t>(siz[1]);
        
        check_arg(m > 0 && n > 0, "m and n must be positive integers.");
        check_arg(I.nelems() == J.nelems(), 
                "I and J must have the same number of elements.");
        
        LMAT_MX_OUT(0, r, marray::numeric_matrix<uint32_t>(m, n), 
                mat_, uint32_t);
        
        dolphin::add_counts(I - T(1), J - T(1), r);
    }
        
}


