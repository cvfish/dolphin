/**********************************************************
 *
 *  _dpaccum.cpp
 *
 *  Internal implementation of dispatched accumulation
 *
 **********************************************************/

#include <light_mat/matlab/matlab_port.h>
#include <dolphin/common/dpaccum.h>
#include <limits>

using namespace lmat;
using namespace lmat::matlab;

const int32_t DPC_SUM = 1;
const int32_t DPC_MAX = 2;
const int32_t DPC_MIN = 3;

LMAT_FP_MEX(dpaccum, 2)
{     
    LMAT_MX(0, code, sca_, int32_t)
    LMAT_MX(2, x, col_, T)
    
    if (nrhs == 4)
    {
        LMAT_MX(1, K, sca_, double)
        LMAT_MX(3, I, col_, int32_t)
        
        const index_t K_ = static_cast<index_t>(K);
        check_arg(K > 0, "K must be a positive integer.");
        
        LMAT_MX_OUT(0, r, marray::numeric_matrix<T>(K_, 1), 
                col_, T);
        
        switch (code)
        {
        case DPC_SUM:
            dolphin::dispatch_sum(x, I, r);
            break;
        case DPC_MAX:
            fill(r, -std::numeric_limits<T>::infinity());
            dolphin::dispatch_max(x, I, r);
            break;
        case DPC_MIN:
            fill(r, std::numeric_limits<T>::infinity());
            dolphin::dispatch_min(x, I, r);
            break;
        }        
    }
    else // nrhs == 5
    {
        LMAT_MX(1, siz, row_, double)
        LMAT_MX(3, I, col_, int32_t)
        LMAT_MX(4, J, col_, int32_t)
                
        check_arg(siz.nelems() == 2, "siz must be a pair of numbers.");
        const index_t m = static_cast<index_t>(siz[0]);
        const index_t n = static_cast<index_t>(siz[1]);
        
        check_arg(m > 0 && n > 0, "m and n must be positive integers.");
        check_arg(I.nelems() == J.nelems(), 
                "I and J must have the same number of elements.");
        
        LMAT_MX_OUT(0, r, marray::numeric_matrix<T>(m, n), 
                mat_, T);
        
        switch (code)
        {
        case DPC_SUM:
            dolphin::dispatch_sum(x, I, J, r);
            break;
        case DPC_MAX:
            fill(r, -std::numeric_limits<T>::infinity());
            dolphin::dispatch_max(x, I, J, r);
            break;
        case DPC_MIN:
            fill(r, std::numeric_limits<T>::infinity());
            dolphin::dispatch_min(x, I, J, r);
            break;
        }  
    }
        
}

