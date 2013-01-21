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

LMAT_FP_MEX(dpaccum, 3)
{     
    LMAT_MX(0, code, sca_, int32_t)
    LMAT_MX(1, dim, sca_, double)
    LMAT_MX(2, K, sca_, double)
    LMAT_MX(3, x, mat_, T)
    LMAT_MX(4, I, col_, int32_t)
    
    const index_t K_ = static_cast<index_t>(K);
    check_arg(K > 0, "K must be a positive integer.");
    
    const index_t m = x.nrows();
    const index_t n = x.ncolumns();    
    
    if (dim == 1) // rows
    {
        LMAT_MX_OUT(0, r, marray::numeric_matrix<T>(K_, n), mat_, T)
        
        switch (code)
        {
        case DPC_SUM:
            dolphin::dispatch_sum_rows(x, I, r);
            break;
        case DPC_MAX:
            fill(r, -std::numeric_limits<T>::infinity());
            dolphin::dispatch_max_rows(x, I, r);
            break;
        case DPC_MIN:
            fill(r, std::numeric_limits<T>::infinity());
            dolphin::dispatch_min_rows(x, I, r);
            break;
        }        
    }
    else // cols
    {
        LMAT_MX_OUT(0, r, marray::numeric_matrix<T>(m, K_), mat_, T)
                        
        switch (code)
        {
        case DPC_SUM:
            dolphin::dispatch_sum_cols(x, I, r);
            break;
        case DPC_MAX:
            fill(r, -std::numeric_limits<T>::infinity());
            dolphin::dispatch_max_cols(x, I, r);
            break;
        case DPC_MIN:
            fill(r, std::numeric_limits<T>::infinity());
            dolphin::dispatch_min_cols(x, I, r);
            break;
        }  
    }
        
}

