/**
 * @file test_common_calc.cpp
 *
 * @brief Test routines in common_calc.h
 *
 * @author Dahua Lin
 */

#include "../test_base.h"
#include <dolphin/common/common_calc.h>

using namespace dolphin;
using namespace dolphin::test;


T_CASE( test_logsumexp )
{
	dense_col<T> x {T(1), T(2), T(3), T(4), T(5)} ;
	T v0 = T(5.45191439593759330506);
	T tol = T(sizeof(T) == 4 ? 1.0e-6 : 1.0e-13);

	exp_terms<T> et;
	et.set_logvalues(x);

	ASSERT_EQ( et.nterms(), 5 );
	ASSERT_EQ( et.max_logvalue(), T(5) );
	ASSERT_APPROX( et.logsum(), v0, tol );
}


template<class A, class B>
void calc_exp(const IRegularMatrix<A, double>& x, IRegularMatrix<B, double>& y)
{
	const index_t n = x.nelems();
	const double *px = x.ptr_data();
	double *py = y.ptr_data();

	typedef lmat::simd_pack<double, lmat::default_simd_kind> pack_t;

	cref_col<double> xv(px, n);
	ref_col<double> yv(py, n);

	lmat::dimension<0> dim(n);

	map(lmat::exp_fun<double>(), lmat::atags::simd<lmat::default_simd_kind>())(4, out_(yv), in_(xv));
/*
	pack_t xp; xp.load_u(px);
	pack_t yp; math::exp(xp);
	yp.store_u(py);
*/
}

T_CASE( test_normalize_exp)
{
	const index_t n = 256;
	dense_col<T> x(n) ;
	for (index_t i = 0; i < n; ++i) x[i] = T(i + 1);

	double se(0);
	for (index_t i = 0; i < n; ++i) se += x[i];

	dense_col<T> p0_(1000);
	ref_col<T> p0(p0_.ptr_data(), n);
	calc_exp(x, p0);

	throw std::runtime_error("runtime error");
}



AUTO_TPACK( test_exp_terms )
{
	ADD_T_CASE( test_logsumexp, float )
	ADD_T_CASE( test_logsumexp, double )
	// ADD_T_CASE( test_normalize_exp, float )
	ADD_T_CASE( test_normalize_exp, double )
}

