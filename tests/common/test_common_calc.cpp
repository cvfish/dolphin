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


T_CASE( test_normalize_exp)
{
	const index_t n = 5;
	dense_col<T> x(n) ;
	for (index_t i = 0; i < n; ++i) x[i] = T(i + 1);

	double se(0);
	for (index_t i = 0; i < n; ++i) se += math::exp(x[i]);

	dense_col<T> p0 = exp(x) / T(se);

	exp_terms<T> et;
	et.set_logvalues(x);

	dense_col<T> p;
	et.normalize_to(p);

	ASSERT_EQ(p.nrows(), 5);
	ASSERT_EQ(p.ncolumns(), 1);

	T tol = T(sizeof(T) == 4 ? 1.0e-6 : 1.0e-13);
	ASSERT_VEC_APPROX(5, p, p0, tol);
}



AUTO_TPACK( test_exp_terms )
{
	ADD_T_CASE( test_logsumexp, float )
	ADD_T_CASE( test_logsumexp, double )
	ADD_T_CASE( test_normalize_exp, float )
	ADD_T_CASE( test_normalize_exp, double )
}

