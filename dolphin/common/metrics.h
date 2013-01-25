/**
 * @file metrics.h
 *
 * @brief Common metrics
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef DOLPHIN_METRICS_H_
#define DOLPHIN_METRICS_H_

#include <dolphin/common/import_lmat.h>
#include <light_mat/linalg/blas_l3.h>
#include <tuple>

#define DOLPHIN_DEF_GENERIC_METRIC(Name, RT) \
	template<typename T> struct Name; \
	template<typename T> \
	struct metric_traits<Name<T> > { \
		typedef T input_type; \
		typedef RT result_type; \
	}; \
	template<typename T> \
	struct Name : public dolphin::IMetric<Name<T> > { \
		typedef T input_type; \
		typedef RT result_type; \
		template<class A, class B> \
		inline RT operator() (const IEWiseMatrix<A, T>& a, const IEWiseMatrix<B, T>& b) const; }; \
	template<typename T> \
	template<class A, class B> \
	inline RT Name<T>::operator() (const IEWiseMatrix<A, T>& a, const IEWiseMatrix<B, T>& b) const


namespace dolphin
{
	/********************************************
	 *
	 *  base class for metrics
	 *
	 ********************************************/

	template<class Metric>
	struct metric_traits;


	template<class Derived>
	class IMetric
	{
	public:
		typedef typename metric_traits<Derived>::result_type result_type;
		typedef typename metric_traits<Derived>::input_type input_type;

		LMAT_CRTP_REF

		template<class A, class B>
		DOLPHIN_ENSURE_INLINE
		result_type operator() (
				const IRegularMatrix<A, result_type>& a,
				const IRegularMatrix<B, result_type>& b) const
		{
			return derived()(a, b);
		}
	};


	template<class Metric, class A, class B, typename TD, class D>
	void colwise(const IMetric<Metric>& metric,
			const IRegularMatrix<A, typename metric_traits<Metric>::input_type>& a,
			const IRegularMatrix<B, typename metric_traits<Metric>::input_type>& b,
			IRegularMatrix<D, TD>& r)
	{
		const index_t na = a.ncolumns();
		const index_t nb = b.ncolumns();

		const Metric& metric_ = metric.derived();
		D& r_ = r.derived();

		if (na == 1)
		{
			for (index_t i = 0; i < nb; ++i)
			{
				r_[i] = static_cast<TD>(metric_(a.column(0), b.column(i)));
			}
		}
		else if (nb == 1)
		{
			for (index_t i = 0; i < na; ++i)
			{
				r_[i] = static_cast<TD>(metric_(a.column(i), b.column(0)));
			}
		}
		else
		{
			LMAT_CHECK_DIMS( na == nb )
			for (index_t i = 0; i < na; ++i)
			{
				r_[i] = static_cast<TD>(metric_(a.column(i), b.column(i)));
			}
		}
	}


	/********************************************
	 *
	 *  pairwise computation
	 *
	 ********************************************/

	template<class Metric, class Arg1, class Arg2>
	class pairwise_metric_expr
	: public IMatrixXpr<pairwise_metric_expr<Metric, Arg1, Arg2>, typename Metric::result_type>
	{
		static const index_t CT_N1 = lmat::meta::ncols<Arg1>::value;
		static const index_t CT_N2 = lmat::meta::ncols<Arg2>::value;
		typedef matrix_shape<CT_N1, CT_N2> shape_type;

	public:
		typedef Metric distance_type;
		typedef Arg1 arg1_type;
		typedef Arg2 arg2_type;

		pairwise_metric_expr(const Metric& metric, const Arg1& a1, const Arg2& a2)
		: m_metric(metric), m_arg1(a1), m_arg2(a2)
		, m_shape(a1.ncolumns(), a2.ncolumns() ){ }

		DOLPHIN_ENSURE_INLINE const Metric& metric() const
		{
			return m_metric;
		}

		DOLPHIN_ENSURE_INLINE const Arg1& arg1() const
		{
			return m_arg1;
		}

		DOLPHIN_ENSURE_INLINE const Arg2& arg2() const
		{
			return m_arg2;
		}

		DOLPHIN_ENSURE_INLINE index_t nrows() const
		{
			return m_shape.nrows();
		}

		DOLPHIN_ENSURE_INLINE index_t ncolumns() const
		{
			return m_shape.ncolumns();
		}

		DOLPHIN_ENSURE_INLINE index_t nelems() const
		{
			return m_shape.nelems();
		}

		DOLPHIN_ENSURE_INLINE shape_type shape() const
		{
			return m_shape;
		}

	private:
		const Metric& m_metric;
		const Arg1& m_arg1;
		const Arg2& m_arg2;
		shape_type m_shape;
	};


	template<class Metric, class Arg1, class Arg2>
	DOLPHIN_ENSURE_INLINE
	inline pairwise_metric_expr<Metric, Arg1, Arg2>
	pairwise(const IMetric<Metric>& metric,
			const IRegularMatrix<Arg1, typename metric_traits<Metric>::input_type>& a1,
			const IRegularMatrix<Arg2, typename metric_traits<Metric>::input_type>& a2)
	{
		return pairwise_metric_expr<Metric, Arg1, Arg2>(metric.derived(), a1.derived(), a2.derived());
	}


	/********************************************
	 *
	 *  metrics between vectors
	 *
	 ********************************************/

	DOLPHIN_DEF_GENERIC_METRIC(euclidean_distance, T)
	{
		return norm(a - b, norms::L2_());
	}

	DOLPHIN_DEF_GENERIC_METRIC(sqeuclidean_distance, T)
	{
		return sqsum(a - b);
	}

	DOLPHIN_DEF_GENERIC_METRIC(cityblock_distance, T)
	{
		return asum(a - b);
	}

	DOLPHIN_DEF_GENERIC_METRIC(chebyshev_distance, T)
	{
		return amax(a - b);
	}


	template<typename T> class minkowski_distance;

	template<typename T>
	struct metric_traits<minkowski_distance<T> >
	{
		typedef T input_type;
		typedef T result_type;
	};

	template<typename T>
	class minkowski_distance : public IMetric<minkowski_distance<T> >
	{
	public:
		DOLPHIN_ENSURE_INLINE
		minkowski_distance(T p)
		: m_p(p), m_inv_p(math::rcp(p)) { }

		DOLPHIN_ENSURE_INLINE
		T p() const { return m_p; }

		template<class A, class B>
		T operator() (const IEWiseMatrix<A, T>& a, const IEWiseMatrix<B, T>& b) const
		{
			return math::pow(sum(pow(abs(a - b), m_p)), m_inv_p);
		}

	private:
		T m_p;
		T m_inv_p;
	};


	DOLPHIN_DEF_GENERIC_METRIC(hamming_distance, uint32_t)
	{
		const index_t n = a.nelems();
		LMAT_CHECK_DIMS( n == b.nelems() )

		auto rd_a = lmat::make_vec_accessor(lmat::scalar_(), in_(a));
		auto rd_b = lmat::make_vec_accessor(lmat::scalar_(), in_(b));

		result_type s(0);
		for (index_t i = 0; i < n; ++i)
			s += result_type(rd_a.scalar(i) != rd_b.scalar(i));
		return s;
	}


	namespace internal
	{
		template<typename T>
		struct cosine_dist_stat
		{
			T xx;
			T xy;
			T yy;

			DOLPHIN_ENSURE_INLINE
			cosine_dist_stat() { }

			DOLPHIN_ENSURE_INLINE
			cosine_dist_stat(const T& x, const T& y)
			: xx(x * x)
			, xy(x * y)
			, yy(y * y)
			{ }

			DOLPHIN_ENSURE_INLINE
			void update(const T& x, const T& y)
			{
				T xx_ = x * x;
				T xy_ = x * y;
				T yy_ = y * y;

				xx += xx_;
				xy += xy_;
				yy += yy_;
			}

			DOLPHIN_ENSURE_INLINE
			void update(const cosine_dist_stat& b)
			{
				xx += b.xx;
				xy += b.xy;
				yy += b.yy;
			}

		};

		template<typename T, typename Kind>
		DOLPHIN_ENSURE_INLINE
		inline cosine_dist_stat<T> reduce_impl(const cosine_dist_stat<lmat::simd_pack<T, Kind> >& s)
		{
			cosine_dist_stat<T> r;
			r.xx = sum(s.xx);
			r.xy = sum(s.xy);
			r.yy = sum(s.yy);
			return r;
		}

		LMAT_DEFINE_AGGREG_SIMD_FOLDKERNEL(cosine_dist_stat, cosine_dist_kernel, 2)
	}

	DOLPHIN_DEF_GENERIC_METRIC(cosine_distance, T)
	{
		const A& a_ = a.derived();
		const B& b_ = b.derived();

		auto r = fold(internal::cosine_dist_kernel<T>())(common_shape(a_, b_), in_(a_), in_(b_));
		return T(1) - ( r.xy / math::sqrt(r.xx * r.yy) );
	}

}


namespace lmat
{
	template<class Dist, class Arg1, class Arg2>
	struct matrix_traits<dolphin::pairwise_metric_expr<Dist, Arg1, Arg2> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = meta::ncols<Arg1>::value;
		static const int ct_num_cols = meta::ncols<Arg2>::value;

		static const bool is_readonly = true;

		typedef matrix_shape<ct_num_rows, ct_num_cols> shape_type;
		typedef typename Dist::result_type value_type;
		typedef typename meta::common_domain<Arg1, Arg2>::type domain;
	};

	LMAT_DEF_SIMD_SUPPORT(dolphin::internal::cosine_dist_kernel)


	// pair metric evaluation

	template<typename Metric, class A, class B, class D>
	void _evaluate(const dolphin::pairwise_metric_expr<Metric, A, B>& expr,
			IRegularMatrix<D, typename dolphin::metric_traits<Metric>::result_type>& dst)
	{
		D& dst_ = dst.derived();
		const index_t m = expr.nrows();
		const index_t n = expr.ncolumns();

		const A& a = expr.arg1();
		const B& b = expr.arg2();

		for (index_t j = 0; j < n; ++j)
		{
			auto bj = b.column(j);

			for (index_t i = 0; i < m; ++i)
			{
				dst_(i, j) = expr.metric()(a.column(i), bj);
			}
		}
	}

	template<class Metric, class A, class B, class D>
	inline void evaluate(const dolphin::pairwise_metric_expr<Metric, A, B>& expr,
			IRegularMatrix<D, typename dolphin::metric_traits<Metric>::result_type>& dst)
	{
		_evaluate(expr, dst);
	}

	template<typename T, class A, class B, class D>
	void evaluate(const dolphin::pairwise_metric_expr<dolphin::sqeuclidean_distance<T>, A, B>& expr,
			IRegularMatrix<D, T>& dst)
	{
		D& dst_ = dst.derived();
		const A& a = expr.arg1();
		const B& b = expr.arg2();

		const index_t m = a.ncolumns();
		const index_t n = b.ncolumns();

		dense_col<T> sa2(m);
		dense_row<T> sb2(n);

		colwise_sqsum(a, sa2);
		colwise_sqsum(b, sb2);

		dst_ = repcol(sa2, n)  + reprow(sb2, m);
		blas::gemm(T(-2), a, b, T(1), dst_, 'T', 'N');
	}

	template<typename T, class A, class B, class D>
	void evaluate(const dolphin::pairwise_metric_expr<dolphin::euclidean_distance<T>, A, B>& expr,
			IRegularMatrix<D, T>& dst)
	{
		D& dst_ = dst.derived();

		dolphin::sqeuclidean_distance<T> sqdist;
		dst_ = dolphin::pairwise(sqdist, expr.arg1(), expr.arg2());
		dst_ = sqrt(dst_);
	}

	template<typename T, class A, class B, class D>
	void evaluate(const dolphin::pairwise_metric_expr<dolphin::cosine_distance<T>, A, B>& expr,
			IRegularMatrix<D, T>& dst)
	{
		D& dst_ = dst.derived();
		const A& a = expr.arg1();
		const B& b = expr.arg2();

		const index_t m = a.ncolumns();
		const index_t n = b.ncolumns();

		dense_col<T> ra(m);
		dense_row<T> rb(n);

		for (index_t i = 0; i < m; ++i) ra[i] = math::rcp(sqsum(a.column(i)));
		for (index_t j = 0; j < n; ++j) rb[j] = math::rcp(sqsum(b.column(j)));

		blas::gemm(a, b, dst_, 'T', 'N');
		dst_ = T(1) - dst_ * sqrt(repcol(ra, n) * reprow(rb, m));
	}
}


#endif
