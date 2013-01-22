/**********************************************************
 *
 *  ktop.cpp
 *
 *  Finding k largest/smallest values
 *
 **********************************************************/

#include <light_mat/matlab/matlab_port.h>
#include <light_mat/mateval/matrix_sort.h>

using namespace lmat;
using namespace lmat::matlab;

template<typename T, typename Ord>
inline void find_ktop(index_t n, index_t k, Ord ord,
        const T *pa, T *ptmp, T *pr)
{
    cref_col<T> a(pa, n);
    ref_col<T> tmp(ptmp, n);
    ref_col<T> r(pr, k);
    
    copy(a, tmp);    
    gsort(tmp, lmat::partial_sort(k), ord);
    copy(tmp(range(0, k)), r);
}


template<typename T, typename Ord>
inline void find_ktop_x(index_t n, index_t k, Ord ord,
        const T *pa, std::pair<T, index_t> *ptmp, T *px, double *pi)
{
    cref_col<T> a(pa, n);
    ref_col<std::pair<T, index_t> > tmp(ptmp, n);

    tmp = gsorted_ex(a, lmat::partial_sort(k), ord);
    
    for (index_t i = 0; i < k; ++i)
    {
        const std::pair<T, index_t>& t = tmp[i];
        px[i] = std::get<0>(t);
        pi[i] = static_cast<double>(std::get<1>(t) + 1);
    }    
}



LMAT_NUM_MEX(ktop, 0)
{
    typedef std::pair<T, index_t> sxv_t;
    
    LMAT_MX_NARGINCHK(2, 2)
    LMAT_MX_NARGOUTCHK(0, 2)
    
    LMAT_MX(0, a, mat_, T)
    LMAT_MX_SCA(1, k, index_t)
    const index_t ak = k >= 0 ? k : -k;
    
    if (is_vector(a))
    {
        const index_t n = a.nelems();
        check_arg(k >= -n && k <= n, "Invalid value of k");
             
        index_t om, on;
        
        if (is_column(a))
        { om = ak; on = 1; }
        else
        { om = 1; on = ak; }
        
        if (nlhs <= 1)
        {
            LMAT_MX_OUT(0, r, marray::numeric_matrix<T>(om, on), col_, T)
            dense_col<T> tmp(n);
            
            if (k >= 0)
            {
                find_ktop(n, ak, asc_(),
                        a.ptr_data(), tmp.ptr_data(), r.ptr_data());
            }
            else
            {
                find_ktop(n, ak, desc_(),
                        a.ptr_data(), tmp.ptr_data(), r.ptr_data());
            }
        }
        else
        {                        
            LMAT_MX_OUT(0, rx, marray::numeric_matrix<T>(om, on), col_, T)
            LMAT_MX_OUT(1, ri, marray::numeric_matrix<double>(om, on), col_, double)
            dense_col<sxv_t> tmp(n);
            
            if (k >= 0)
            {
                find_ktop_x(n, ak, asc_(),
                        a.ptr_data(), tmp.ptr_data(), 
                        rx.ptr_data(), ri.ptr_data());
            }
            else
            {
                find_ktop_x(n, ak, desc_(),
                        a.ptr_data(), tmp.ptr_data(), 
                        rx.ptr_data(), ri.ptr_data());
            }
        }

    }
    else
    {
        const index_t m = a.nrows();
        const index_t n = a.ncolumns();
        check_arg(k >= -m && k <= m, "Invalid value of k");
                
        if (nlhs <= 1)
        {
            LMAT_MX_OUT(0, r, marray::numeric_matrix<T>(ak, n), mat_, T)
            dense_col<T> tmp(m);
            
            if (k >= 0)
            {
                for (index_t j = 0; j < n; ++j)
                {
                    find_ktop(m, ak, asc_(),
                            a.ptr_col(j), tmp.ptr_data(), r.ptr_col(j));
                }
            }
            else
            {
                mexPrintf("desc\n");
                for (index_t j = 0; j < n; ++j)
                {
                    find_ktop(m, ak, desc_(),
                            a.ptr_col(j), tmp.ptr_data(), r.ptr_col(j));
                }
            }
        }
        else
        {
            LMAT_MX_OUT(0, rx, marray::numeric_matrix<T>(ak, n), mat_, T)
            LMAT_MX_OUT(1, ri, marray::numeric_matrix<double>(ak, n), mat_, double)
            dense_col<sxv_t> tmp(m);
            
            if (k >= 0)
            {
                for (index_t j = 0; j < n; ++j)
                {
                    find_ktop_x(m, ak, asc_(),
                            a.ptr_col(j), tmp.ptr_data(), 
                            rx.ptr_col(j), ri.ptr_col(j));
                }
            }
            else
            {
                mexPrintf("desc\n");
                for (index_t j = 0; j < n; ++j)
                {
                    find_ktop_x(m, ak, desc_(),
                            a.ptr_col(j), tmp.ptr_data(), 
                            rx.ptr_col(j), ri.ptr_col(j));
                }
            }
        }
    }
    
}

