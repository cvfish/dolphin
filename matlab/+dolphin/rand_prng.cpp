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

rand_stream_handle _default_rstream()
{
    static random::default_rand_stream rstream;
    return reinterpret_cast<rand_stream_handle>(&rstream);
}


LMAT_SIMPLE_MEX( rand_prng )
{
    LMAT_MX_NARGINCHK(0, 0)
    LMAT_MX_NARGOUTCHK(0, 1)
    
    plhs[0] = marray::from_scalar(_default_rstream());
}

