/**********************************************************
 *
 *  symeig.cpp
 *
 *  symmetric eigenvalue analysis
 *
 **********************************************************/

#include <light_mat/matlab/matlab_port.h>
#include <light_mat/linalg/lapack_syev.h>
#include <algorithm>

using namespace lmat;
using namespace lmat::matlab;
using std::swap;

inline char check_op(char c)
{
    if (c == 'R' || c == 'r') return 'R';
    else if (c == 'D' || c == 'd') return 'D';
    else if (c == 'N' || c == 'n') return 'N';
    else
        throw invalid_argument("Invalid value for c");
}


template<typename T>
inline void flip(ref_col<T>& evs)
{
    index_t l = 0;
    index_t r = evs.nelems() - 1;
    while (l < r)
    {
        swap(evs[l], evs[r]);
        ++l; --r;
    }
}


template<typename T>
inline void flip(ref_col<T>& evs, ref_matrix<T>& U)
{
    flip(evs);
    
    index_t n = U.ncolumns();
    index_t l = 0;
    index_t r = n - 1;
    
    dense_col<T> tmp(n);
    
    while (l < r)
    {
        auto cl = U.column(l);
        auto cr = U.column(r);
        
        copy(cl, tmp);
        copy(cr, cl);
        copy(tmp, cr);
        
        ++l;
        --r;
    }
}



LMAT_FP_MEX(symeig, 0)
{
    LMAT_MX_NARGINCHK(1, 2)
    LMAT_MX_NARGOUTCHK(0, 2)
    
    LMAT_MX(0, a, mat_, T)
    
    char op = 'R';
    if (nrhs == 2)
    {
        LMAT_MX_SCA(1, c, char)
        op = check_op(c);
    }
    
    
    check_arg(is_square(a), "a must be a square matrix.");
    const index_t n = a.nrows();
    
    LMAT_MX_OUT(0, evs, marray::numeric_matrix<T>(n, 1), col_, T)
    
    if (op == 'N')
    {
        if (nlhs <= 1)
        {
            lapack::syev(a, evs);
            flip(evs);
        }
        else
        {
            LMAT_MX_OUT(1, U, marray::numeric_matrix<T>(n, n), mat_, T)
            lapack::syev(a, evs, U);
            flip(evs, U);
        }
    }
    else if (op == 'R')
    {                
        if (nlhs <= 1)
        {
            lapack::syevr(a, evs);
            flip(evs);
        }
        else
        {
            LMAT_MX_OUT(1, U, marray::numeric_matrix<T>(n, n), mat_, T)
            lapack::syevr(a, evs, U);
            flip(evs, U);
        }
    }
    else if (op == 'D')
    {                
        if (nlhs <= 1)
        {
            lapack::syevd(a, evs);
            flip(evs);
        }
        else
        {
            LMAT_MX_OUT(1, U, marray::numeric_matrix<T>(n, n), mat_, T)
            lapack::syevd(a, evs, U);
            flip(evs, U);
        }
    }
}



