/**********************************************************
 *
 *  rand_prng.cpp
 *
 *  get a handle to Dolphin's random number generator
 *
 **********************************************************/

#include "rand_base.h"

using namespace lmat;
using namespace lmat::matlab;

random::default_rand_stream* _ptr_default_rstream()
{
    static random::default_rand_stream rstream;
    return &rstream;
}


LMAT_SIMPLE_MEX( rand_prng )
{
    LMAT_MX_NARGINCHK(0, 1)
    LMAT_MX_NARGOUTCHK(0, 1)
    
    random::default_rand_stream *p = _ptr_default_rstream();
    
    if (nrhs >= 1)
    {
        LMAT_MX(0, seed, sca_, double)
        p->set_seed(static_cast<uint32_t>(seed));
    }
    
    if (nlhs >= 1)
        plhs[0] = marray::from_scalar(
                reinterpret_cast<rand_stream_handle>(p));
}

