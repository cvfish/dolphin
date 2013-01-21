/**********************************************************
 *
 * pd_logdet.cpp
 *
 * Log-determinant of positive definite matrix
 *
 **********************************************************/

#include <light_mat/matlab/matlab_port.h>
#include <light_mat/linalg/lapack_chol.h>

using namespace lmat;
using namespace lmat::matlab;

LMAT_FP_MEX(pd_logdet, 0)
{
    LMAT_MX_NARGINCHK(1, 1)
    LMAT_MX_NARGOUTCHK(0, 1)
    
    LMAT_MX(0, a, mat_, T)
    
    T v = pdlogdet(a);

    plhs[0] = marray::from_scalar(v);
    
}