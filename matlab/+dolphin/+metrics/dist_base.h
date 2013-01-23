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

template<typename T, class Dist>
void distance_port(
        mxArray *plhs[],
        const Dist& dist,
        const cref_matrix<T>& a, 
        const cref_matrix<T>& b)
{
    typedef typename Dist::result_type RT;
    
    if (is_vector(a))
    {
        if (is_vector(b))
        {
            check_arg(a.nelems() == b.nelems(), 
                    "Inconsistent dimensions of a and b.");
            
            RT r = dist(as_col(a), as_col(b));
            plhs[0] = marray::from_scalar(r);
        }
        else
        {
            check_arg(a.nelems() == b.nrows(), 
                    "Inconsistent dimensions of a and b.");
            
            const index_t n = b.ncolumns();
            LMAT_MX_OUT(0, r, marray::numeric_matrix<RT>(1, n), row_, RT);
            dolphin::colwise(dist, as_col(a), b, r);
        }
    }
    else
    {
        if (is_vector(b))
        {            
            check_arg(a.nrows() == b.nelems(), 
                    "Inconsistent dimensions of a and b.");
            
            const index_t n = a.ncolumns();                        
            LMAT_MX_OUT(0, r, marray::numeric_matrix<RT>(1, n), row_, RT);
            dolphin::colwise(dist, a, as_col(b), r);
        }
        else
        {
            check_arg(have_same_shape(a, b),
                    "Inconsistent dimensions of a and b.");
                                    
            const index_t n = a.ncolumns();                        
            LMAT_MX_OUT(0, r, marray::numeric_matrix<RT>(1, n), row_, RT);
            dolphin::colwise(dist, a, b, r);
        }
    }
}
