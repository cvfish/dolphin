/**********************************************************
 *
 *  The basis for implementing matlab distance functions
 *
 **********************************************************/

#include <light_mat/matlab/matlab_port.h>
#include <light_mat/matrix/matrix_asvec.h>
#include <dolphin/common/metrics.h>

using namespace lmat;
using namespace lmat::matlab;


inline void _dist_check_dim(index_t rdim, index_t da, index_t db)
{
    if (rdim)
    {
        check_arg(da == rdim && db == rdim,
                "Inconsistent dimensions of a, b, and w");
    }
    else
    {
        check_arg(da == db,
                "Inconsistent dimensions of a and b.");
    }
}


template<typename T, class Dist>
void distance_port(
        mxArray *plhs[],
        const Dist& dist,
        const cref_matrix<T>& a, 
        const cref_matrix<T>& b, 
        index_t require_dim = 0)
{
    typedef typename Dist::result_type RT;
    
    if (is_vector(a))
    {
        if (is_vector(b))
        {
            _dist_check_dim(require_dim, a.nelems(), b.nelems());
            
            RT r = dist(as_col(a), as_col(b));
            plhs[0] = marray::from_scalar(r);
        }
        else
        {                        
            _dist_check_dim(require_dim, a.nelems(), b.nrows());
            
            const index_t n = b.ncolumns();
            LMAT_MX_OUT(0, r, marray::numeric_matrix<RT>(1, n), row_, RT);
            dolphin::colwise(dist, as_col(a), b, r);
        }
    }
    else
    {
        if (is_vector(b))
        {            
            _dist_check_dim(require_dim, a.nrows(), b.nelems());
            
            const index_t n = a.ncolumns();                        
            LMAT_MX_OUT(0, r, marray::numeric_matrix<RT>(1, n), row_, RT);
            dolphin::colwise(dist, a, as_col(b), r);
        }
        else
        {
            _dist_check_dim(require_dim, a.nrows(), b.nrows());                                            
            const index_t n = a.ncolumns();
            check_arg( n == b.ncolumns(), 
                    "a and b have different number of columns");
            
            LMAT_MX_OUT(0, r, marray::numeric_matrix<RT>(1, n), row_, RT);
            dolphin::colwise(dist, a, b, r);
        }
    }
}


template<typename T, class Dist>
void pw_distance_port(        
        mxArray *plhs[],
        const Dist& dist,
        const cref_matrix<T>& a, 
        const cref_matrix<T>& b, 
        index_t require_dim = 0)
{
    typedef typename Dist::result_type RT;
    
    if (is_empty(b))
    {
        _dist_check_dim(require_dim, a.nrows(), a.nrows());
        
        const index_t n = a.ncolumns();
        
        LMAT_MX_OUT(0, r, marray::numeric_matrix<RT>(n, n), mat_, RT);
        r = dolphin::pairwise(dist, a);
    }
    else
    {                
        _dist_check_dim(require_dim, a.nrows(), b.nrows());
        
        const index_t m = a.ncolumns();
        const index_t n = b.ncolumns();
        
        LMAT_MX_OUT(0, r, marray::numeric_matrix<RT>(m, n), mat_, RT);
        r = dolphin::pairwise(dist, a, b);
    }
}






