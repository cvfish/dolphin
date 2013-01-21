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

	dense_col<T> p0(n);
	for (index_t i = 0; i < n; ++i) p0[i] = math::exp(x[i]) / T(se);

	exp_terms<T> et;
	et.set_logvalues(x);

	dense_col<T> p;
	et.normalize_to(p);

	ASSERT_EQ(p.nrows(), n);
	ASSERT_EQ(p.ncolumns(), 1);

	T tol = T(sizeof(T) == 4 ? 1.0e-6 : 1.0e-13);
	ASSERT_VEC_APPROX(n, p, p0, tol);
}

T_CASE( test_full_entropy )
{
	const index_t n = 10;
	dense_col<T> p(n);
	fill_randr(p, T(-0.5), T(1.0));

	double v0(0);
	for (index_t i = 0; i < n; ++i)
	{
		if (p[i] > 0) v0 -= p[i] * math::log(p[i]);
	}

	T v = entropy(p);

	T tol = T(sizeof(T) == 4 ? 1.0e-6 : 1.0e-13);
	ASSERT_APPROX(v, v0, tol);
}

T_CASE( test_colwise_entropy )
{
	const index_t m = 10;
	const index_t n = 12;
	dense_matrix<T> p(m, n);
	fill_randr(p, T(-0.2), T(1.0));

	dense_row<T> r(n);
	colwise_entropy(p, r);

	dense_row<T> r0(n);
	for (index_t j = 0; j < n; ++j)
		r0[j] = entropy(p.column(j));

	T tol = T(sizeof(T) == 4 ? 1.0e-6 : 1.0e-13);
	ASSERT_VEC_APPROX(m, r, r0, tol);
}


T_CASE( test_rowwise_entropy )
{
	const index_t m = 10;
	const index_t n = 12;
	dense_matrix<T> p(m, n);
	fill_randr(p, T(-0.2), T(1.0));

	dense_row<T> r(m);
	rowwise_entropy(p, r);

	dense_row<T> r0(m);
	for (index_t i = 0; i < m; ++i)
		r0[i] = entropy(p.row(i));

	T tol = T(sizeof(T) == 4 ? 1.0e-6 : 1.0e-13);
	ASSERT_VEC_APPROX(m, r, r0, tol);
}


AUTO_TPACK( test_exp_terms )
{
	ADD_T_CASE( test_logsumexp, float )
	ADD_T_CASE( test_logsumexp, double )
	ADD_T_CASE( test_normalize_exp, float )
	ADD_T_CASE( test_normalize_exp, double )
}


AUTO_TPACK( test_entropy )
{
	ADD_T_CASE( test_full_entropy, float )
	ADD_T_CASE( test_full_entropy, double )
	ADD_T_CASE( test_colwise_entropy, float )
	ADD_T_CASE( test_colwise_entropy, double )
	ADD_T_CASE( test_rowwise_entropy, float )
	ADD_T_CASE( test_rowwise_entropy, double )
}



