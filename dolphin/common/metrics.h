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

#define DOLPHIN_DEF_GENERIC_METRIC_TRAITS(Name, RT, IsPosDef, IsSym) \
	template<typename T> \
	struct metric_traits<Name<T> > { \
		typedef T input_type; \
		typedef RT result_type; \
		static const bool is_positive_definite = IsPosDef; \
		static const bool is_symmetric = IsSym; \
	};

#define DOLPHIN_DEF_GENERIC_WMETRIC_TRAITS(Name, RT, IsPosDef, IsSym) \
	template<typename T, typename W> \
	struct metric_traits<Name<T, W> > { \
		typedef T input_type; \
		typedef RT result_type; \
		static const bool is_positive_definite = IsPosDef; \
		static const bool is_symmetric = IsSym; \
	};

#define DOLPHIN_DEF_GENERIC_METRIC(Name, RT, IsPosDef, IsSym) \
	template<typename T> class Name; \
	DOLPHIN_DEF_GENERIC_METRIC_TRAITS(Name, RT, IsPosDef, IsSym) \
	template<typename T> \
	class Name : public dolphin::IMetric<Name<T> > { \
	public: \
		typedef T input_type; \
		typedef RT result_type; \
		template<class A, class B> \
		inline RT operator() (const IEWiseMatrix<A, T>& a, const IEWiseMatrix<B, T>& b) const; }; \
	template<typename T> \
	template<class A, class B> \
	inline RT Name<T>::operator() (const IEWiseMatrix<A, T>& a, const IEWiseMatrix<B, T>& b) const

#define DOLPHIN_DEF_GENERIC_WMETRIC(Name, RT, IsPosDef, IsSym) \
	template<typename T, typename W> class Name; \
	DOLPHIN_DEF_GENERIC_WMETRIC_TRAITS(Name, RT, IsPosDef, IsSym) \
	template<typename T, typename W> \
	class Name : public dolphin::IMetric<Name<T, W> > { \
	public: \
		typedef T input_type; \
		typedef RT weight_type; \
		typedef RT result_type; \
		DOLPHIN_ENSURE_INLINE \
		Name(const W& w) : m_weights(w) { } \
		DOLPHIN_ENSURE_INLINE \
		const W& weights() const { return m_weights; } \
		template<class A, class B> \
		inline RT operator() (const IEWiseMatrix<A, T>& a, const IEWiseMatrix<B, T>& b) const; \
	private: \
		const W& m_weights; \
	}; \
	template<typename T, typename W> \
	template<class A, class B> \
	inline RT Name<T, W>::operator() (const IEWiseMatrix<A, T>& a, const IEWiseMatrix<B, T>& b) const




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
	: public lmat::matrix_xpr_base<pairwise_metric_expr<Metric, Arg1, Arg2> >
	{
		typedef lmat::matrix_xpr_base<pairwise_metric_expr<Metric, Arg1, Arg2> > base_t;

	public:
		typedef Metric metric_type;
		typedef Arg1 arg1_type;
		typedef Arg2 arg2_type;

		pairwise_metric_expr(const Metric& metric, const Arg1& a1, const Arg2& a2)
		: base_t(a1.ncolumns(), a2.ncolumns())
		, m_metric(metric), m_arg1(a1), m_arg2(a2) { }

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

	private:
		const Metric& m_metric;
		const Arg1& m_arg1;
		const Arg2& m_arg2;
	};

	template<class Metric, class Arg>
	class self_pairwise_metric_expr
	: public lmat::matrix_xpr_base<self_pairwise_metric_expr<Metric, Arg> >
	{
		typedef lmat::matrix_xpr_base<self_pairwise_metric_expr<Metric, Arg> > base_t;

	public:
		typedef Metric metric_type;
		typedef Arg arg_type;

		self_pairwise_metric_expr(const Metric& metric, const Arg& a)
		: base_t(a.ncolumns(), a.ncolumns())
		, m_metric(metric), m_arg(a) { }

		DOLPHIN_ENSURE_INLINE const Metric& metric() const
		{
			return m_metric;
		}

		DOLPHIN_ENSURE_INLINE const Arg& arg() const
		{
			return m_arg;
		}

	private:
		const Metric& m_metric;
		const Arg& m_arg;
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


	template<class Metric, class Arg>
	DOLPHIN_ENSURE_INLINE
	inline self_pairwise_metric_expr<Metric, Arg>
	pairwise(const IMetric<Metric>& metric,
			const IRegularMatrix<Arg, typename metric_traits<Metric>::input_type>& a)
	{
		return self_pairwise_metric_expr<Metric, Arg>(metric.derived(), a.derived());
	}


	/********************************************
	 *
	 *  metrics between vectors
	 *
	 ********************************************/

	// Euclidean

	DOLPHIN_DEF_GENERIC_METRIC(euclidean_distance, T, true, true)
	{
		return norm(a - b, norms::L2_());
	}

	DOLPHIN_DEF_GENERIC_WMETRIC(weuclidean_distance, T, true, true)
	{
		return math::sqrt(sum(weights() * sqr(a - b)));
	}

	template<typename T, typename W>
	DOLPHIN_ENSURE_INLINE
	inline weuclidean_distance<T, W> weighted_euclidean(const IRegularMatrix<W, T>& weights)
	{
		LMAT_CHECK_DIMS( is_column(weights) );
		return weuclidean_distance<T, W>(weights.derived());
	}

	// Squared Euclidean

	DOLPHIN_DEF_GENERIC_METRIC(sqeuclidean_distance, T, true, true)
	{
		return sqsum(a - b);
	}

	DOLPHIN_DEF_GENERIC_WMETRIC(wsqeuclidean_distance, T, true, true)
	{
		return sum(weights() * sqr(a - b));
	}

	template<typename T, typename W>
	DOLPHIN_ENSURE_INLINE
	inline wsqeuclidean_distance<T, W> weighted_sqeuclidean(const IRegularMatrix<W, T>& weights)
	{
		LMAT_CHECK_DIMS( is_column(weights) );
		return wsqeuclidean_distance<T, W>(weights.derived());
	}

	// City block

	DOLPHIN_DEF_GENERIC_METRIC(cityblock_distance, T, true, true)
	{
		return asum(a - b);
	}

	DOLPHIN_DEF_GENERIC_WMETRIC(wcityblock_distance, T, true, true)
	{
		return sum(weights() * abs(a - b));
	}

	template<typename T, typename W>
	DOLPHIN_ENSURE_INLINE
	inline wcityblock_distance<T, W> weighted_cityblock(const IRegularMatrix<W, T>& weights)
	{
		LMAT_CHECK_DIMS( is_column(weights) );
		return wcityblock_distance<T, W>(weights.derived());
	}

	// Chebyshev

	DOLPHIN_DEF_GENERIC_METRIC(chebyshev_distance, T, true, true)
	{
		return amax(a - b);
	}


	// Minkowski

	template<typename T> class minkowski_distance;
	DOLPHIN_DEF_GENERIC_METRIC_TRAITS(minkowski_distance, T, true, true)

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

	template<typename T, typename W> class wminkowski_distance;
	DOLPHIN_DEF_GENERIC_WMETRIC_TRAITS(wminkowski_distance, T, true, true)

	template<typename T, typename W>
	class wminkowski_distance : public IMetric<wminkowski_distance<T, W> >
	{
	public:
		DOLPHIN_ENSURE_INLINE
		wminkowski_distance(T p, const W& weights)
		: m_p(p), m_inv_p(math::rcp(p)), m_weights(weights) { }

		DOLPHIN_ENSURE_INLINE
		T p() const { return m_p; }

		DOLPHIN_ENSURE_INLINE
		const W& weights() const
		{
			return m_weights;
		}

		template<class A, class B>
		T operator() (const IEWiseMatrix<A, T>& a, const IEWiseMatrix<B, T>& b) const
		{
			return math::pow(sum(weights() * pow(abs(a - b), m_p)), m_inv_p);
		}

	private:
		T m_p;
		T m_inv_p;
		const W& m_weights;
	};


	template<typename T, typename W>
	DOLPHIN_ENSURE_INLINE
	inline wminkowski_distance<T, W> weighted_minkowski(const T& p, const IRegularMatrix<W, T>& weights)
	{
		LMAT_CHECK_DIMS( is_column(weights) );
		return wminkowski_distance<T, W>(p, weights.derived());
	}


	// Hamming

	DOLPHIN_DEF_GENERIC_METRIC(hamming_distance, uint32_t, true, true)
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

	DOLPHIN_DEF_GENERIC_WMETRIC(whamming_distance, typename lmat::matrix_traits<W>::value_type, true, true)
	{
		const index_t n = a.nelems();
		LMAT_CHECK_DIMS( n == b.nelems() )

		auto rd_a = lmat::make_vec_accessor(lmat::scalar_(), in_(a));
		auto rd_b = lmat::make_vec_accessor(lmat::scalar_(), in_(b));

		const W& w_ = weights();

		result_type s(0);
		for (index_t i = 0; i < n; ++i)
			s += result_type(rd_a.scalar(i) != rd_b.scalar(i)) * w_[i];
		return s;
	}

	template<typename T, typename W, typename TW>
	DOLPHIN_ENSURE_INLINE
	inline whamming_distance<T, W> weighted_hamming(type_<T>, const IRegularMatrix<W, TW>& weights)
	{
		LMAT_CHECK_DIMS( is_column(weights) );
		return whamming_distance<T, W>(weights.derived());
	}


	// cosine distance

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

	DOLPHIN_DEF_GENERIC_METRIC(cosine_distance, T, true, true)
	{
		const A& a_ = a.derived();
		const B& b_ = b.derived();

		auto r = fold(internal::cosine_dist_kernel<T>())(common_shape(a_, b_), in_(a_), in_(b_));
		return T(1) - ( r.xy / math::sqrt(r.xx * r.yy) );
	}

}


namespace lmat
{
	// matrix traits for pairwise metrics

	template<class Dist, class Arg1, class Arg2>
	struct matrix_traits<dolphin::pairwise_metric_expr<Dist, Arg1, Arg2> >
	: public lmat::matrix_xpr_traits_base<
	  typename Dist::result_type,
	  meta::ncols<Arg1>::value,
	  meta::ncols<Arg2>::value,
	  typename meta::common_domain<Arg1, Arg2>::type> { };

	template<class Dist, class Arg>
	struct matrix_traits<dolphin::self_pairwise_metric_expr<Dist, Arg> >
	: public lmat::matrix_xpr_traits_base<
	  typename Dist::result_type,
	  meta::ncols<Arg>::value,
	  meta::ncols<Arg>::value,
	  typename meta::domain_of<Arg>::type> { };


	// introduce SIMD support to cosine_dist_kernel

	LMAT_DEF_SIMD_SUPPORT(dolphin::internal::cosine_dist_kernel)


	/********************************************
	 *
	 *  generic pairwise evaluation
	 *
	 ********************************************/

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

	template<typename Metric, class A, class D>
	void _evaluate(const dolphin::self_pairwise_metric_expr<Metric, A>& expr,
			IRegularMatrix<D, typename dolphin::metric_traits<Metric>::result_type>& dst)
	{
		D& dst_ = dst.derived();
		const index_t n = expr.ncolumns();

		const A& a = expr.arg();

		const bool is_pos_def = dolphin::metric_traits<Metric>::is_positive_definite;
		const bool is_symmetric = dolphin::metric_traits<Metric>::is_symmetric;
		typedef typename dolphin::metric_traits<Metric>::result_type RT;

		if (is_symmetric)
		{
			if (is_pos_def)
			{
				for (index_t j = 0; j < n; ++j)
				{
					auto aj = a.column(j);

					for (index_t i = 0; i < j; ++i)
						dst_(i, j) = dst_(j, i);

					dst_(j, j) = RT(0);

					for (index_t i = j+1; i < n; ++i)
						dst_(i, j) = expr.metric()(a.column(i), aj);
				}
			}
			else
			{
				for (index_t j = 0; j < n; ++j)
				{
					auto aj = a.column(j);

					for (index_t i = 0; i < j; ++i)
						dst_(i, j) = dst_(j, i);

					for (index_t i = j; i < n; ++i)
						dst_(i, j) = expr.metric()(a.column(i), aj);
				}
			}
		}
		else
		{
			if (is_pos_def)
			{
				for (index_t j = 0; j < n; ++j)
				{
					auto aj = a.column(j);

					for (index_t i = 0; i < j; ++i)
						dst_(i, j) = expr.metric()(a.column(i), aj);

					dst_(j, j) = RT(0);

					for (index_t i = j+1; i < n; ++i)
						dst_(i, j) = expr.metric()(a.column(i), aj);
				}
			}
			else
			{
				for (index_t j = 0; j < n; ++j)
				{
					auto aj = a.column(j);

					for (index_t i = 0; i < n; ++i)
						dst_(i, j) = expr.metric()(a.column(i), aj);
				}
			}
		}
	}


	template<class Metric, class A, class B, class D>
	inline void evaluate(const dolphin::pairwise_metric_expr<Metric, A, B>& expr,
			IRegularMatrix<D, typename dolphin::metric_traits<Metric>::result_type>& dst)
	{
		_evaluate(expr, dst);
	}

	template<class Metric, class A, class D>
	inline void evaluate(const dolphin::self_pairwise_metric_expr<Metric, A>& expr,
			IRegularMatrix<D, typename dolphin::metric_traits<Metric>::result_type>& dst)
	{
		_evaluate(expr, dst);
	}


	/********************************************
	 *
	 *  specialized pairwise evaluation
	 *
	 ********************************************/

	// sqeuclidean_distance

	template<typename T, class D>
	DOLPHIN_ENSURE_INLINE
	inline void _postprocess_posdef_metrics(IRegularMatrix<D, T>& dst, bool selfpw)
	{
		D& dst_ = dst.derived();

		dst_ = max(dst_, T(0));

		if (selfpw)
		{
			auto dvs = dst_.diag();
			dvs << T(0);
		}
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

		_postprocess_posdef_metrics(dst_, false);
	}

	template<typename T, class A, class D>
	void evaluate(const dolphin::self_pairwise_metric_expr<dolphin::sqeuclidean_distance<T>, A>& expr,
			IRegularMatrix<D, T>& dst)
	{
		D& dst_ = dst.derived();
		const A& a = expr.arg();

		const index_t n = a.ncolumns();

		dense_col<T> sa2(n);
		colwise_sqsum(a, sa2);

		dst_ = repcol(sa2, n)  + reprow(as_row(sa2), n);
		blas::gemm(T(-2), a, a, T(1), dst_, 'T', 'N');

		_postprocess_posdef_metrics(dst_, true);
	}

	template<typename T, class W, class A, class B, class D>
	void evaluate(const dolphin::pairwise_metric_expr<dolphin::wsqeuclidean_distance<T, W>, A, B>& expr,
			IRegularMatrix<D, T>& dst)
	{
		D& dst_ = dst.derived();
		const A& a = expr.arg1();
		const B& b = expr.arg2();
		const W& w = expr.metric().weights();

		const index_t m = a.ncolumns();
		const index_t n = b.ncolumns();

		dense_col<T> sa2(m);
		dense_row<T> sb2(n);

		colwise_sum(sqr(a) * repcol(w, m), sa2);
		colwise_sum(sqr(b) * repcol(w, n), sb2);

		dst_ = repcol(sa2, n)  + reprow(sb2, m);

		dense_matrix<T> wa = a * repcol(w, m);
		blas::gemm(T(-2), wa, b, T(1), dst_, 'T', 'N');

		_postprocess_posdef_metrics(dst_, false);
	}

	template<typename T, class W, class A, class D>
	void evaluate(const dolphin::self_pairwise_metric_expr<dolphin::wsqeuclidean_distance<T, W>, A>& expr,
			IRegularMatrix<D, T>& dst)
	{
		D& dst_ = dst.derived();
		const A& a = expr.arg();
		const W& w = expr.metric().weights();

		const index_t n = a.ncolumns();

		dense_col<T> sa2(n);
		colwise_sum(sqr(a) * repcol(w, n), sa2);

		dst_ = repcol(sa2, n)  + reprow(as_row(sa2), n);

		dense_matrix<T> wa = a * repcol(w, n);
		blas::gemm(T(-2), wa, a, T(1), dst_, 'T', 'N');

		_postprocess_posdef_metrics(dst_, true);
	}

	// euclidean_distance

	template<typename T, class A, class B, class D>
	void evaluate(const dolphin::pairwise_metric_expr<dolphin::euclidean_distance<T>, A, B>& expr,
			IRegularMatrix<D, T>& dst)
	{
		D& dst_ = dst.derived();

		dolphin::sqeuclidean_distance<T> sqdist;
		dst_ = dolphin::pairwise(sqdist, expr.arg1(), expr.arg2());
		dst_ = sqrt(dst_);
	}

	template<typename T, class A, class D>
	void evaluate(const dolphin::self_pairwise_metric_expr<dolphin::euclidean_distance<T>, A>& expr,
			IRegularMatrix<D, T>& dst)
	{
		D& dst_ = dst.derived();

		dolphin::sqeuclidean_distance<T> sqdist;
		dst_ = dolphin::pairwise(sqdist, expr.arg());
		dst_ = sqrt(dst_);
	}

	template<typename T, class W, class A, class B, class D>
	void evaluate(const dolphin::pairwise_metric_expr<dolphin::weuclidean_distance<T, W>, A, B>& expr,
			IRegularMatrix<D, T>& dst)
	{
		D& dst_ = dst.derived();

		dolphin::wsqeuclidean_distance<T, W> sqdist(expr.metric().weights());
		dst_ = dolphin::pairwise(sqdist, expr.arg1(), expr.arg2());
		dst_ = sqrt(dst_);
	}

	template<typename T, class W, class A, class D>
	void evaluate(const dolphin::self_pairwise_metric_expr<dolphin::weuclidean_distance<T, W>, A>& expr,
			IRegularMatrix<D, T>& dst)
	{
		D& dst_ = dst.derived();

		dolphin::wsqeuclidean_distance<T, W> sqdist(expr.metric().weights());
		dst_ = dolphin::pairwise(sqdist, expr.arg());
		dst_ = sqrt(dst_);
	}


	// cosine_distance

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

		_postprocess_posdef_metrics(dst_, false);
	}

	template<typename T, class A, class D>
	void evaluate(const dolphin::self_pairwise_metric_expr<dolphin::cosine_distance<T>, A>& expr,
			IRegularMatrix<D, T>& dst)
	{
		D& dst_ = dst.derived();
		const A& a = expr.arg();

		const index_t n = a.ncolumns();

		dense_col<T> ra(n);

		for (index_t i = 0; i < n; ++i) ra[i] = math::rcp(sqsum(a.column(i)));

		blas::gemm(a, a, dst_, 'T', 'N');
		dst_ = T(1) - dst_ * sqrt(repcol(ra, n) * reprow(as_row(ra), n));

		_postprocess_posdef_metrics(dst_, true);
	}
}


#endif
