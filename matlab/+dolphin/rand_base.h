// declaration of external functions to facilitate random simulation

#include <light_mat/matlab/matlab_port.h>
#include <light_mat/random/distr_fwd.h>

template<size_t PtrSize>
struct _rand_stream_handle_helper;

template<>
struct _rand_stream_handle_helper<4>
{
    typedef lmat::uint32_t type;
};

template<>
struct _rand_stream_handle_helper<8>
{
    typedef lmat::uint64_t type;
};

typedef typename _rand_stream_handle_helper<
        sizeof(lmat::random::default_rand_stream*)
        >::type rand_stream_handle;

inline lmat::random::default_rand_stream* get_default_rstream()
{
    mxArray *mx = 0;
    mexCallMATLAB(1, &mx, 0, 0, "dolphin.rand_prng");
    
    rand_stream_handle h = *((const rand_stream_handle*)mxGetData(mx));
    return reinterpret_cast<lmat::random::default_rand_stream*>(h);
}

