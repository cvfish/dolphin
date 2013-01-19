/**
 * @file test_dpaccum.cpp
 *
 * Unit testing of dispatched accumulation
 * 
 * @author Dahua Lin 
 */

#include "../test_base.h"
#include <dolphin/common/dpaccum.h>

using namespace dolphin;
using namespace dolphin::test;

SIMPLE_CASE( test_add_counts )
{
	const index_t K = 12;
	const index_t len = 200;
	dense_col<uint32_t> cnts(K, zero());
	dense_col<index_t> I(len);

	fill_randi(I, (index_t)0, K+2);

	dense_col<uint32_t> c0(K, zero());
	for (index_t i = 0; i < len; ++i)
	{
		index_t k = I[i];
		if (k >= 0 && k < K) c0[k]++;
	}

	add_counts(I, cnts);
	ASSERT_VEC_EQ(K, cnts, c0);
}


SIMPLE_CASE( test_dispatch_sum_1d )
{
	const index_t K = 12;
	const index_t len = 200;

	dense_col<index_t> I(len);
	dense_col<double> v(len);
	dense_col<double> a(K, zero());
	dense_col<double> a0(K, zero());

	fill_randi(I, (index_t)0, K+2);
	fill_randr(v, 0., 1.);

	for (index_t i = 0; i < len; ++i)
	{
		index_t k = I[i];
		if (k >= 0 && k < K) a0[k] += v[i];
	}

	dispatch_sum(v, I, a);

	ASSERT_VEC_APPROX(K, a, a0, 1.0e-15);
}

SIMPLE_CASE( test_dispatch_max_1d )
{
	const index_t K = 12;
	const index_t len = 200;

	dense_col<index_t> I(len);
	dense_col<double> v(len);
	dense_col<double> a(K, zero());
	dense_col<double> a0(K, zero());

	fill_randi(I, (index_t)0, K+2);
	fill_randr(v, 0., 1.);

	for (index_t i = 0; i < len; ++i)
	{
		index_t k = I[i];
		if (k >= 0 && k < K)
		{
			if (v[i] > a0[k]) a0[k] = v[i];
		}
	}

	dispatch_max(v, I, a);

	ASSERT_VEC_EQ(K, a, a0);
}


SIMPLE_CASE( test_dispatch_min_1d )
{
	const index_t K = 12;
	const index_t len = 200;

	dense_col<index_t> I(len);
	dense_col<double> v(len);
	dense_col<double> a(K, zero());
	dense_col<double> a0(K, zero());

	fill_randi(I, (index_t)0, K+2);
	fill_randr(v, -2., 0.);

	for (index_t i = 0; i < len; ++i)
	{
		index_t k = I[i];
		if (k >= 0 && k < K)
		{
			if (v[i] < a0[k]) a0[k] = v[i];
		}
	}

	dispatch_min(v, I, a);

	ASSERT_VEC_EQ(K, a, a0);
}


SIMPLE_CASE( test_dispatch_sum_2d )
{
	const index_t M = 3;
	const index_t N = 4;
	const index_t len = 200;

	dense_col<index_t> I(len);
	dense_col<index_t> J(len);
	dense_col<double> v(len);
	dense_matrix<double> a(M, N, zero());
	dense_matrix<double> a0(M, N, zero());

	fill_randi(I, (index_t)0, M+1);
	fill_randi(J, (index_t)0, N+1);
	fill_randr(v, 0., 1.);

	for (index_t i = 0; i < len; ++i)
	{
		index_t ci = I[i];
		index_t cj = J[i];
		if (ci >= 0 && ci < M && cj >= 0 && cj < N)
		{
			a0(ci, cj) += v[i];
		}
	}

	dispatch_sum(v, I, J, a);

	ASSERT_MAT_APPROX(M, N, a, a0, 1.0e-15);
}

SIMPLE_CASE( test_dispatch_max_2d )
{
	const index_t M = 3;
	const index_t N = 4;
	const index_t len = 200;

	dense_col<index_t> I(len);
	dense_col<index_t> J(len);
	dense_col<double> v(len);
	dense_matrix<double> a(M, N, zero());
	dense_matrix<double> a0(M, N, zero());

	fill_randi(I, (index_t)0, M+1);
	fill_randi(J, (index_t)0, N+1);
	fill_randr(v, 0., 1.);

	for (index_t i = 0; i < len; ++i)
	{
		index_t ci = I[i];
		index_t cj = J[i];
		if (ci >= 0 && ci < M && cj >= 0 && cj < N)
		{
			if (v[i] > a0(ci, cj)) a0(ci, cj) = v[i];
		}
	}

	dispatch_max(v, I, J, a);

	ASSERT_MAT_APPROX(M, N, a, a0, 1.0e-15);
}

SIMPLE_CASE( test_dispatch_min_2d )
{
	const index_t M = 3;
	const index_t N = 4;
	const index_t len = 200;

	dense_col<index_t> I(len);
	dense_col<index_t> J(len);
	dense_col<double> v(len);
	dense_matrix<double> a(M, N, zero());
	dense_matrix<double> a0(M, N, zero());

	fill_randi(I, (index_t)0, M+1);
	fill_randi(J, (index_t)0, N+1);
	fill_randr(v, -2., 0.);

	for (index_t i = 0; i < len; ++i)
	{
		index_t ci = I[i];
		index_t cj = J[i];
		if (ci >= 0 && ci < M && cj >= 0 && cj < N)
		{
			if (v[i] < a0(ci, cj)) a0(ci, cj) = v[i];
		}
	}

	dispatch_min(v, I, J, a);

	ASSERT_MAT_APPROX(M, N, a, a0, 1.0e-15);
}


SIMPLE_CASE( test_dispatch_sum_cols )
{
	const index_t m = 15;
	const index_t n = 60;
	const index_t K = 5;

	dense_col<index_t> J(n);
	dense_matrix<double> v(m, n);
	dense_matrix<double> a(m, K, zero());
	dense_matrix<double> a0(m, K, zero());

	fill_randi(J, (index_t)0, K+1);
	fill_randr(v, 0., 1.);

	for (index_t j = 0; j < n; ++j)
	{
		index_t cj = J[j];

		if (cj >= 0 && cj < K)
		{
			auto rj = a0.column(cj);
			rj += v.column(j);
		}
	}

	dispatch_sum_cols(v, J, a);

	ASSERT_MAT_APPROX(m, K, a, a0, 1.0e-15);
}


SIMPLE_CASE( test_dispatch_max_cols )
{
	const index_t m = 15;
	const index_t n = 60;
	const index_t K = 5;

	dense_col<index_t> J(n);
	dense_matrix<double> v(m, n);
	dense_matrix<double> a(m, K, zero());
	dense_matrix<double> a0(m, K, zero());

	fill_randi(J, (index_t)0, K+1);
	fill_randr(v, 0., 1.);

	for (index_t j = 0; j < n; ++j)
	{
		index_t cj = J[j];

		if (cj >= 0 && cj < K)
		{
			auto rj = a0.column(cj);
			rj = max(rj, v.column(j));
		}
	}

	dispatch_max_cols(v, J, a);

	ASSERT_MAT_APPROX(m, K, a, a0, 1.0e-15);
}


SIMPLE_CASE( test_dispatch_min_cols )
{
	const index_t m = 15;
	const index_t n = 60;
	const index_t K = 5;

	dense_col<index_t> J(n);
	dense_matrix<double> v(m, n);
	dense_matrix<double> a(m, K, zero());
	dense_matrix<double> a0(m, K, zero());

	fill_randi(J, (index_t)0, K+1);
	fill_randr(v, -2., 0.);

	for (index_t j = 0; j < n; ++j)
	{
		index_t cj = J[j];

		if (cj >= 0 && cj < K)
		{
			auto rj = a0.column(cj);
			rj = min(rj, v.column(j));
		}
	}

	dispatch_min_cols(v, J, a);

	ASSERT_MAT_APPROX(m, K, a, a0, 1.0e-15);
}


SIMPLE_CASE( test_dispatch_sum_rows )
{
	const index_t m = 60;
	const index_t n = 15;
	const index_t K = 5;

	dense_col<index_t> I(m);
	dense_matrix<double> v(m, n);
	dense_matrix<double> a(K, n, zero());
	dense_matrix<double> a0(K, n, zero());

	fill_randi(I, (index_t)0, K+1);
	fill_randr(v, 0., 1.);

	for (index_t i = 0; i < m; ++i)
	{
		index_t ci = I[i];

		if (ci >= 0 && ci < K)
		{
			auto ri = a0.row(ci);
			ri += v.row(i);
		}
	}

	dispatch_sum_rows(v, I, a);
	ASSERT_MAT_APPROX(K, n, a, a0, 1.0e-15);
}


SIMPLE_CASE( test_dispatch_max_rows )
{
	const index_t m = 60;
	const index_t n = 15;
	const index_t K = 5;

	dense_col<index_t> I(m);
	dense_matrix<double> v(m, n);
	dense_matrix<double> a(K, n, zero());
	dense_matrix<double> a0(K, n, zero());

	fill_randi(I, (index_t)0, K+1);
	fill_randr(v, 0., 1.);

	for (index_t i = 0; i < m; ++i)
	{
		index_t ci = I[i];

		if (ci >= 0 && ci < K)
		{
			auto ri = a0.row(ci);
			ri = max(ri, v.row(i));
		}
	}

	dispatch_max_rows(v, I, a);
	ASSERT_MAT_APPROX(K, n, a, a0, 1.0e-15);
}


SIMPLE_CASE( test_dispatch_min_rows )
{
	const index_t m = 60;
	const index_t n = 15;
	const index_t K = 5;

	dense_col<index_t> I(m);
	dense_matrix<double> v(m, n);
	dense_matrix<double> a(K, n, zero());
	dense_matrix<double> a0(K, n, zero());

	fill_randi(I, (index_t)0, K+1);
	fill_randr(v, -2., 0.);

	for (index_t i = 0; i < m; ++i)
	{
		index_t ci = I[i];

		if (ci >= 0 && ci < K)
		{
			auto ri = a0.row(ci);
			ri = min(ri, v.row(i));
		}
	}

	dispatch_min_rows(v, I, a);
	ASSERT_MAT_APPROX(K, n, a, a0, 1.0e-15);
}



AUTO_TPACK( test_counts )
{
	ADD_SIMPLE_CASE( test_add_counts )
}

AUTO_TPACK( test_dpaccum_1d )
{
	ADD_SIMPLE_CASE( test_dispatch_sum_1d )
	ADD_SIMPLE_CASE( test_dispatch_max_1d )
	ADD_SIMPLE_CASE( test_dispatch_min_1d )
}

AUTO_TPACK( test_dpaccum_2d )
{
	ADD_SIMPLE_CASE( test_dispatch_sum_2d )
	ADD_SIMPLE_CASE( test_dispatch_max_2d )
	ADD_SIMPLE_CASE( test_dispatch_min_2d )
}

AUTO_TPACK( test_dpaccum_cols )
{
	ADD_SIMPLE_CASE( test_dispatch_sum_cols )
	ADD_SIMPLE_CASE( test_dispatch_max_cols )
	ADD_SIMPLE_CASE( test_dispatch_min_cols )
}


AUTO_TPACK( test_dpaccum_rows )
{
	ADD_SIMPLE_CASE( test_dispatch_sum_rows )
	ADD_SIMPLE_CASE( test_dispatch_max_rows )
	ADD_SIMPLE_CASE( test_dispatch_min_rows )
}



