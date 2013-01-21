/**
 * @file common_calc.h
 *
 * @brief Common calculation routines
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef DOLPHIN_COMMON_CALC_H_
#define DOLPHIN_COMMON_CALC_H_

#include <dolphin/common/import_lmat.h>

namespace dolphin
{

	template<typename T>
	class exp_terms
	{
		static_assert(std::is_floating_point<T>::value,
				"T must be floating-point types.");

	public:
		DOLPHIN_ENSURE_INLINE
		explicit exp_terms()
		: m_tmp()
		, m_max_logval(0)
		, m_shift(0)
		, m_tmp_sum(0) { }

		template<class A>
		void set_logvalues(const IRegularMatrix<A, T>& a)
		{
			m_max_logval = maximum(a);
			m_shift = m_max_logval;
			m_tmp = exp(a - m_shift);
			m_tmp_sum = sum(m_tmp);
		}

		DOLPHIN_ENSURE_INLINE
		index_t nterms() const
		{
			return m_tmp.nelems();
		}

		DOLPHIN_ENSURE_INLINE
		T max_logvalue() const
		{
			return m_max_logval;
		}

		DOLPHIN_ENSURE_INLINE
		T logsum() const
		{
			return m_shift + math::log(m_tmp_sum);
		}

		template<class P>
		void normalize_to(IRegularMatrix<P, T>& p)
		{
			T c = math::rcp(m_tmp_sum);
			p.derived() = m_tmp * c;
		}

	private:
		dense_matrix<T> m_tmp;  // stores exp(x - shift)
		T m_max_logval;
		T m_shift;
		T m_tmp_sum;
	};


	template<typename T, class P>
	inline T entropy(const IEWiseMatrix<P, T>& p)
	{
		return - sum(xlogx(p));
	}

	template<typename T, class P, class R>
	inline void colwise_entropy(const IEWiseMatrix<P, T>& p, IRegularMatrix<R, T>& r)
	{
		colwise_sum(xlogx(p), r);
		r.derived() = -r;
	}

	template<typename T, class P, class R>
	inline void rowwise_entropy(const IEWiseMatrix<P, T>& p, IRegularMatrix<R, T>& r)
	{
		rowwise_sum(xlogx(p), r);
		r.derived() = -r;
	}


}

#endif
