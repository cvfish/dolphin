/**
 * @file test_metrics.cpp
 *
 * @brief Unit testing of common metrics
 *
 * @author Dahua Lin
 */


#include "../test_base.h"
#include <dolphin/common/metrics.h>
#include <functional>

using namespace dolphin;
using namespace dolphin::test;

typedef dense_matrix<double> mat_t;
typedef cref_matrix<double, 0, 1> col_t;

const index_t vdim = 13;
const index_t M = 7;
const index_t N = 8;


/************************************************
 *
 *  safe implementation as reference
 *
 ************************************************/

template<class A, class B, typename T, typename Fun>
dense_matrix<typename Fun::result_type>
my_pairwise(const IRegularMatrix<A, T>& a, const IRegularMatrix<B, T>& b, const Fun& fun)
{
	const index_t m = a.ncolumns();
	const index_t n = b.ncolumns();

	dense_matrix<typename Fun::result_type> dists(m, n);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			dists(i, j) = fun(a.column(i), b.column(j));
		}
	}

	return dists;
}

template<class A, class B, class W, typename T, typename Fun>
dense_matrix<typename Fun::result_type>
my_wpairwise(const IRegularMatrix<A, T>& a, const IRegularMatrix<B, T>& b, const IRegularMatrix<W, T>& w, const Fun& fun)
{
	const index_t m = a.ncolumns();
	const index_t n = b.ncolumns();

	dense_matrix<typename Fun::result_type> dists(m, n);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			dists(i, j) = fun(a.column(i), b.column(j), w.derived());
		}
	}

	return dists;
}


#define DEF_MY_DIST(Name) \
		struct Name { \
			typedef double T; \
			typedef double result_type; \
			double operator()(const col_t& a, const col_t& b) const; }; \
		double Name::operator() (const col_t& a, const col_t& b) const

#define DEF_MY_WDIST(Name) \
		struct Name { \
			typedef double T; \
			typedef double result_type; \
			double operator()(const col_t& a, const col_t& b, const dense_col<double>& w) const; }; \
		double Name::operator() (const col_t& a, const col_t& b, const dense_col<double>& w) const

DEF_MY_DIST( my_euclidean_distance )
{
	T d(0);
	for (index_t i = 0; i < a.nelems(); ++i)
	{
		T v = a[i] - b[i];
		d += v * v;
	}
	return math::sqrt(d);
}

DEF_MY_WDIST( my_weighted_euclidean )
{
	T d(0);
	for (index_t i = 0; i < a.nelems(); ++i)
	{
		T v = a[i] - b[i];
		d += w[i] * v * v;
	}
	return math::sqrt(d);
}


DEF_MY_DIST( my_sqeuclidean_distance )
{
	T d(0);
	for (index_t i = 0; i < a.nelems(); ++i)
	{
		T v = a[i] - b[i];
		d += v * v;
	}
	return d;
}

DEF_MY_WDIST( my_weighted_sqeuclidean )
{
	T d(0);
	for (index_t i = 0; i < a.nelems(); ++i)
	{
		T v = a[i] - b[i];
		d += w[i] * v * v;
	}
	return d;
}


DEF_MY_DIST( my_cityblock_distance )
{
	T d(0);
	for (index_t i = 0; i < a.nelems(); ++i)
	{
		d += math::abs(a[i] - b[i]);
	}
	return d;
}

DEF_MY_WDIST( my_weighted_cityblock )
{
	T d(0);
	for (index_t i = 0; i < a.nelems(); ++i)
	{
		d += w[i] * math::abs(a[i] - b[i]);
	}
	return d;
}


DEF_MY_DIST( my_chebyshev_distance )
{
	T d(0);
	for (index_t i = 0; i < a.nelems(); ++i)
	{
		T v = math::abs(a[i] - b[i]);
		if (v > d) d = v;
	}
	return d;
}


DEF_MY_DIST( my_minkowski_distance )
{
	T d(0);
	T e = T(3.2);
	for (index_t i = 0; i < a.nelems(); ++i)
	{
		T v = math::abs(a[i] - b[i]);
		d += math::pow(v, e);
	}
	return math::pow(d, math::rcp(e));
}

DEF_MY_WDIST( my_weighted_minkowski )
{
	T d(0);
	T e = T(3.2);
	for (index_t i = 0; i < a.nelems(); ++i)
	{
		T v = math::abs(a[i] - b[i]);
		d += w[i] * math::pow(v, e);
	}
	return math::pow(d, math::rcp(e));
}


DEF_MY_DIST( my_cosine_distance )
{
	T xx(0), yy(0), xy(0);
	for (index_t i = 0; i < a.nelems(); ++i)
	{
		xx += a[i] * a[i];
		xy += a[i] * b[i];
		yy += b[i] * b[i];
	}
	return T(1) - xy / math::sqrt(xx * yy);
}


#define DEF_DIST_TEST_(Name, Construct) \
		SIMPLE_CASE( test_##Name ) { \
			mat_t a(vdim, M); \
			mat_t b(vdim, N); \
			fill_randr(a, -1.0, 1.0); \
			fill_randr(b, -1.0, 1.0); \
			Construct; \
			mat_t D0 = my_pairwise(a, b, my_##Name()); \
			mat_t D1 = my_pairwise(a, b, dist); \
			mat_t D2 = pairwise(dist, a, b); \
			ASSERT_EQ( D1.nrows(), M ); \
			ASSERT_EQ( D1.ncolumns(), N ); \
			ASSERT_EQ( D2.nrows(), M ); \
			ASSERT_EQ( D2.ncolumns(), N ); \
			double tol = 1.0e-13; \
			ASSERT_MAT_APPROX(M, N, D1, D0, tol); \
			ASSERT_MAT_APPROX(M, N, D2, D0, tol); \
			mat_t S0 = my_pairwise(a, a, my_##Name()); \
			mat_t S2 = pairwise(dist, a); \
			ASSERT_EQ( S2.nrows(), M ); \
			ASSERT_EQ( S2.ncolumns(), M ); \
			ASSERT_MAT_APPROX(M, M, S2, S0, tol); }

#define DEF_WDIST_TEST_(Name, Construct) \
		SIMPLE_CASE( test_##Name ) { \
			mat_t a(vdim, M); \
			mat_t b(vdim, N); \
			dense_col<double> w(vdim); \
			fill_randr(a, -1.0, 1.0); \
			fill_randr(b, -1.0, 1.0); \
			fill_randr(w, 0.0, 2.0); \
			Construct; \
			mat_t D0 = my_wpairwise(a, b, w, my_##Name()); \
			mat_t D1 = my_pairwise(a, b, dist); \
			mat_t D2 = pairwise(dist, a, b); \
			ASSERT_EQ( D1.nrows(), M ); \
			ASSERT_EQ( D1.ncolumns(), N ); \
			ASSERT_EQ( D2.nrows(), M ); \
			ASSERT_EQ( D2.ncolumns(), N ); \
			double tol = 1.0e-13; \
			ASSERT_MAT_APPROX(M, N, D1, D0, tol); \
			ASSERT_MAT_APPROX(M, N, D2, D0, tol); \
			mat_t S0 = my_wpairwise(a, a, w, my_##Name()); \
			mat_t S2 = pairwise(dist, a); \
			ASSERT_EQ( S2.nrows(), M ); \
			ASSERT_EQ( S2.ncolumns(), M ); \
			ASSERT_MAT_APPROX(M, M, S2, S0, tol); }

#define DEF_DIST_TEST(Name) DEF_DIST_TEST_( Name, Name<double> dist )
#define DEF_WDIST_TEST(Name) DEF_WDIST_TEST_( Name, auto dist = Name(w) )

DEF_DIST_TEST( euclidean_distance )
DEF_DIST_TEST( sqeuclidean_distance )
DEF_DIST_TEST( cityblock_distance )
DEF_DIST_TEST( chebyshev_distance )
DEF_DIST_TEST_( minkowski_distance, minkowski_distance<double> dist(3.2) )

DEF_DIST_TEST( cosine_distance )

SIMPLE_CASE( test_hamming_distance )
{
	mat_t a(vdim, M);
	mat_t b(vdim, N);

	fill_randi(a, 1., 3.);
	fill_randi(b, 1., 3.);

	dense_matrix<uint32_t> D0(M, N);

	for (index_t j = 0; j < N; ++j)
	{
		for (index_t i = 0; i < M; ++i)
		{
			uint32_t s0(0);
			for (index_t k = 0; k < vdim; ++k)
				if (a(k, i) != b(k, j)) ++ s0;

			D0(i, j) = s0;
		}
	}

	dense_matrix<uint32_t> D1 = my_pairwise(a, b, hamming_distance<double>());

	ASSERT_MAT_EQ( M, N, D1, D0 );

	dense_matrix<uint32_t> D2 = pairwise(hamming_distance<double>(), a, b);

	ASSERT_EQ( D2.nrows(), M );
	ASSERT_EQ( D2.ncolumns(), N );
	ASSERT_MAT_EQ( M, N, D2, D0 );
}

// weighted distances

DEF_WDIST_TEST( weighted_euclidean )
DEF_WDIST_TEST( weighted_sqeuclidean )
DEF_WDIST_TEST( weighted_cityblock )

inline wminkowski_distance<double, dense_matrix<double, 0, 1> > wminkowski_(const dense_col<double>& w)
{
	return weighted_minkowski(3.2, w);
}

DEF_WDIST_TEST_( weighted_minkowski, auto dist = wminkowski_(w) )


SIMPLE_CASE( test_weighted_hamming )
{
	dense_matrix<int32_t> a(vdim, M);
	dense_matrix<int32_t> b(vdim, N);
	dense_col<double> w(vdim);

	fill_randi(a, 1, 3);
	fill_randi(b, 1, 3);
	fill_randr(w, 0.0, 2.0);

	dense_matrix<double> D0(M, N);

	for (index_t j = 0; j < N; ++j)
	{
		for (index_t i = 0; i < M; ++i)
		{
			double s0(0);
			for (index_t k = 0; k < vdim; ++k)
				if (a(k, i) != b(k, j)) s0 += w[k];

			D0(i, j) = s0;
		}
	}

	dense_matrix<double> D1 = my_pairwise(a, b, weighted_hamming(type_<int32_t>(), w));

	ASSERT_MAT_APPROX( M, N, D1, D0, 1.0e-12 );

	dense_matrix<double> D2 = pairwise(weighted_hamming(type_<int32_t>(), w), a, b);

	ASSERT_EQ( D2.nrows(), M );
	ASSERT_EQ( D2.ncolumns(), N );
	ASSERT_MAT_APPROX( M, N, D2, D0, 1.0e-12 );
}



// colwise evaluation

SIMPLE_CASE( colwise_metric_00 )
{
	mat_t a(vdim, M);
	mat_t b(vdim, M);

	fill_randr(a, 0.0, 1.0);
	fill_randr(b, 0.0, 1.0);

	sqeuclidean_distance<double> dist;

	dense_matrix<double> D0(1, M);
	for (index_t i = 0; i < M; ++i)
	{
		D0[i] = dist(a.column(i), b.column(i));
	}

	dense_matrix<double> D1(1, M);
	colwise(dist, a, b, D1);

	double tol = 1.0e-14;
	ASSERT_VEC_APPROX(M, D0, D1, tol);
}

SIMPLE_CASE( colwise_metric_01 )
{
	mat_t a(vdim, M);
	mat_t b(vdim, 1);

	fill_randr(a, 0.0, 1.0);
	fill_randr(b, 0.0, 1.0);

	sqeuclidean_distance<double> dist;

	dense_matrix<double> D0(1, M);
	for (index_t i = 0; i < M; ++i)
	{
		D0[i] = dist(a.column(i), b);
	}

	dense_matrix<double> D1(1, M);
	colwise(dist, a, b, D1);

	double tol = 1.0e-14;
	ASSERT_VEC_APPROX(M, D0, D1, tol);
}

SIMPLE_CASE( colwise_metric_10 )
{
	mat_t a(vdim, 1);
	mat_t b(vdim, M);

	fill_randr(a, 0.0, 1.0);
	fill_randr(b, 0.0, 1.0);

	sqeuclidean_distance<double> dist;

	dense_matrix<double> D0(1, M);
	for (index_t i = 0; i < M; ++i)
	{
		D0[i] = dist(a, b.column(i));
	}

	dense_matrix<double> D1(1, M);
	colwise(dist, a, b, D1);

	double tol = 1.0e-14;
	ASSERT_VEC_APPROX(M, D0, D1, tol);
}



AUTO_TPACK( basic_dists )
{
	ADD_SIMPLE_CASE( test_euclidean_distance )
	ADD_SIMPLE_CASE( test_sqeuclidean_distance )
	ADD_SIMPLE_CASE( test_cityblock_distance )
	ADD_SIMPLE_CASE( test_chebyshev_distance )
	ADD_SIMPLE_CASE( test_minkowski_distance )

	ADD_SIMPLE_CASE( test_cosine_distance )
	ADD_SIMPLE_CASE( test_hamming_distance )
}

AUTO_TPACK( weighted_dists )
{
	ADD_SIMPLE_CASE( test_weighted_euclidean )
	ADD_SIMPLE_CASE( test_weighted_sqeuclidean )
	ADD_SIMPLE_CASE( test_weighted_cityblock )
	ADD_SIMPLE_CASE( test_weighted_minkowski )

	ADD_SIMPLE_CASE( test_weighted_hamming )
}

AUTO_TPACK( colwise_dists )
{
	ADD_SIMPLE_CASE( colwise_metric_00 )
	ADD_SIMPLE_CASE( colwise_metric_01 )
	ADD_SIMPLE_CASE( colwise_metric_10 )
}



