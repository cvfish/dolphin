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
    LMAT_MX_SCA(0, code, int32_t)
    LMAT_MX_SCA(1, dim, index_t)
    LMAT_MX_SCA(2, K, index_t)
    LMAT_MX(3, x, mat_, T)
    LMAT_MX(4, I, col_, int32_t)
    
    check_arg(K > 0, "K must be a positive integer.");
    
    const index_t m = x.nrows();
    const index_t n = x.ncolumns();    
    
    if (dim == 1) // rows
    {
        LMAT_MX_OUT(0, r, marray::numeric_matrix<T>(K, n), mat_, T)
        
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
        LMAT_MX_OUT(0, r, marray::numeric_matrix<T>(m, K), mat_, T)
                        
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

