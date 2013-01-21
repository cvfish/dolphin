/**
 * @file dpaccum.h
 *
 * Functions for dispatched accumulation
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_DPACCUM_H_
#define LIGHTMAT_DPACCUM_H_

#include <dolphin/common/import_lmat.h>
#include <light_mat/mateval/mat_reduce.h>

namespace dolphin
{
	template<typename TI, class Indices, typename TC, class Counts>
	inline typename std::enable_if<
		lmat::supports_linear_macc<Indices>::value &&
		supports_linear_index<Counts>::value,
	void>::type
	_add_counts(
			const IEWiseMatrix<Indices, TI>& I,
			IRegularMatrix<Counts, TC>& counts)
	{
		Counts& cnts = counts.derived();

		const index_t n = I.nelems();
		const index_t K = counts.nelems();
		auto rd = lmat::make_vec_accessor(lmat::atags::scalar(), in_(I.derived()));

		for (index_t i = 0; i < n; ++i)
		{
			index_t k = static_cast<index_t>(rd.scalar(i));

			if (k >= 0 && k < K)
			{
				++ cnts[k];
			}
		}
	}

	template<typename TI, class Iinds, typename TJ, class Jinds, typename TC, class Counts>
	inline typename std::enable_if<
		lmat::supports_linear_macc<Iinds>::value &&
		lmat::supports_linear_macc<Jinds>::value,
	void>::type
	_add_counts(
			const IEWiseMatrix<Iinds, TI>& I,
			const IEWiseMatrix<Jinds, TJ>& J,
			IRegularMatrix<Counts, TC>& counts)
	{
		Counts& cnts = counts.derived();

		const index_t n = I.nelems();
		const index_t M = counts.nrows();
		const index_t N = counts.ncolumns();
		auto rd_i = lmat::make_vec_accessor(lmat::atags::scalar(), in_(I.derived()));
		auto rd_j = lmat::make_vec_accessor(lmat::atags::scalar(), in_(J.derived()));

		for (index_t i = 0; i < n; ++i)
		{
			index_t ci = static_cast<index_t>(rd_i.scalar(i));
			index_t cj = static_cast<index_t>(rd_j.scalar(i));

			if (ci >= 0 && ci < M && cj >= 0 && cj < N)
			{
				++ cnts(ci, cj);
			}
		}
	}


	template<typename TI, class Indices, typename TC, class Counts>
	inline void add_counts(
			const IEWiseMatrix<Indices, TI>& I,
			IRegularMatrix<Counts, TC>& counts)
	{
		_add_counts(I, counts);
	}


	template<typename TI, class Iinds, typename TJ, class Jinds, typename TC, class Counts>
	inline void add_counts(
			const IEWiseMatrix<Iinds, TI>& I,
			const IEWiseMatrix<Jinds, TJ>& J,
			IRegularMatrix<Counts, TC>& counts)
	{
		_add_counts(I, J, counts);
	}


	template<typename TI, class ISubs, typename T, class Values, class Result, class Folder>
	inline typename std::enable_if<
		lmat::supports_linear_macc<ISubs>::value &&
		lmat::supports_linear_macc<Values>::value &&
		supports_linear_index<Result>::value,
	void>::type
	dispatch_accum(
			const IEWiseMatrix<Values, T>& values,
			const IEWiseMatrix<ISubs, TI>& I,
			IRegularMatrix<Result, T>& result,
			const Folder& folder)
	{
		const Values& v = values.derived();
		Result& r = result.derived();

		const index_t n = v.nelems();
		const index_t K = r.nelems();

		check_arg(I.nelems() == n, "The sizes of I and values are inconsistent.");

		auto rd_l = lmat::make_vec_accessor(lmat::atags::scalar(), in_(I.derived()));
		auto rd_v = lmat::make_vec_accessor(lmat::atags::scalar(), in_(v));

		for (index_t i = 0; i < n; ++i)
		{
			index_t k = static_cast<index_t>(rd_l.scalar(i));

			if (k >= 0 && k < K)
			{
				folder.fold(r[k], rd_v.scalar(i));
			}
		}
	}


	template<typename TI, class ISubs, typename T, class Values, class Result>
	inline void dispatch_sum(
			const IEWiseMatrix<Values, T>& values,
			const IEWiseMatrix<ISubs, TI>& I,
			IRegularMatrix<Result, T>& result)
	{
		dispatch_accum(values, I, result, lmat::sum_folder<T>());
	}

	template<typename TI, class ISubs, typename T, class Values, class Result>
	inline void dispatch_max(
			const IEWiseMatrix<Values, T>& values,
			const IEWiseMatrix<ISubs, TI>& I,
			IRegularMatrix<Result, T>& result)
	{
		dispatch_accum(values, I, result, lmat::maximum_folder<T>());
	}

	template<typename TI, class ISubs, typename T, class Values, class Result>
	inline void dispatch_min(
			const IEWiseMatrix<Values, T>& values,
			const IEWiseMatrix<ISubs, TI>& I,
			IRegularMatrix<Result, T>& result)
	{
		dispatch_accum(values, I, result, lmat::minimum_folder<T>());
	}


	template<typename TI, class ISubs, class JSubs,
		typename T, class Values, class Result, class Folder>
	inline typename std::enable_if<
		lmat::supports_linear_macc<ISubs>::value &&
		lmat::supports_linear_macc<JSubs>::value &&
		lmat::supports_linear_macc<Values>::value &&
		supports_linear_index<Result>::value,
	void>::type
	dispatch_accum(
			const IEWiseMatrix<Values, T>& values,
			const IEWiseMatrix<ISubs, TI>& I,
			const IEWiseMatrix<JSubs, TI>& J,
			IRegularMatrix<Result, T>& result,
			const Folder& folder)
	{
		const Values& v = values.derived();
		Result& r = result.derived();

		const index_t n = v.nelems();
		const index_t M = r.nrows();
		const index_t N = r.ncolumns();

		check_arg(I.nelems() == n && J.nelems() == n,
				"The sizes of I, J, and values are inconsistent.");

		auto rd_i = lmat::make_vec_accessor(lmat::atags::scalar(), in_(I.derived()));
		auto rd_j = lmat::make_vec_accessor(lmat::atags::scalar(), in_(J.derived()));
		auto rd_v = lmat::make_vec_accessor(lmat::atags::scalar(), in_(v));

		for (index_t i = 0; i < n; ++i)
		{
			index_t ci = static_cast<index_t>(rd_i.scalar(i));
			index_t cj = static_cast<index_t>(rd_j.scalar(i));

			if (ci >= 0 && ci < M && cj >= 0 && cj < N)
			{
				folder.fold(r(ci, cj), rd_v.scalar(i));
			}
		}
	}

	template<typename TI, class ISubs, class JSubs, typename T, class Values, class Result>
	inline void dispatch_sum(
			const IEWiseMatrix<Values, T>& values,
			const IEWiseMatrix<ISubs, TI>& I,
			const IEWiseMatrix<JSubs, TI>& J,
			IRegularMatrix<Result, T>& result)
	{
		dispatch_accum(values, I, J, result, lmat::sum_folder<T>());
	}

	template<typename TI, class ISubs, class JSubs, typename T, class Values, class Result>
	inline void dispatch_max(
			const IEWiseMatrix<Values, T>& values,
			const IEWiseMatrix<ISubs, TI>& I,
			const IEWiseMatrix<JSubs, TI>& J,
			IRegularMatrix<Result, T>& result)
	{
		dispatch_accum(values, I, J, result, lmat::maximum_folder<T>());
	}

	template<typename TI, class ISubs, class JSubs, typename T, class Values, class Result>
	inline void dispatch_min(
			const IEWiseMatrix<Values, T>& values,
			const IEWiseMatrix<ISubs, TI>& I,
			const IEWiseMatrix<JSubs, TI>& J,
			IRegularMatrix<Result, T>& result)
	{
		dispatch_accum(values, I, J, result, lmat::minimum_folder<T>());
	}


	template<typename TI, class JSubs, typename T, class Values, class Result, class Folder>
	inline typename std::enable_if<
		lmat::supports_linear_macc<JSubs>::value &&
		is_percol_contiguous<Values>::value &&
		is_percol_contiguous<Result>::value,
	void>::type
	dispatch_accum_cols(
			const IEWiseMatrix<Values, T>& values,
			const IEWiseMatrix<JSubs, TI>& J,
			IRegularMatrix<Result, T>& result,
			const Folder& folder)
	{
		const Values& v = values.derived();
		Result& r = result.derived();

		const index_t m = v.nrows();
		const index_t n = v.ncolumns();
		const index_t K = r.ncolumns();

		check_arg( J.nelems() == n, "The sizes of J and values are inconsistent" );
		check_arg( r.nrows() == m, "The numbers of rows in values and result are inconsistent." );

		auto rd_j = lmat::make_vec_accessor(lmat::atags::scalar(), in_(J.derived()));

		lmat::_parfold_kernel<Folder> fker(folder);

		typedef lmat::default_simd_kind skind;

		typedef typename std::conditional<
				lmat::is_simdizable<Folder, skind>::value &&
				lmat::supports_simd<Values, skind, false>::value &&
				lmat::supports_simd<Values, skind, false>::value,
				lmat::atags::simd<lmat::default_simd_kind>,
				lmat::atags::scalar
		>::type atag;

		auto vker = ewise(fker, atag());

		auto rd = lmat::make_multicol_accessor(atag(), in_(v));
		auto wt = lmat::make_multicol_accessor(atag(), in_out_(r));

		for (index_t j = 0; j < n; ++j)
		{
			index_t cj = static_cast<index_t>(rd_j.scalar(j));
			if (cj >= 0 && cj < K)
			{
				vker.apply(m, wt.col(cj), rd.col(j));
			}
		}
	}


	template<typename TI, class JSubs, typename T, class Values, class Result>
	inline void dispatch_sum_cols(
			const IEWiseMatrix<Values, T>& values,
			const IEWiseMatrix<JSubs, TI>& J,
			IRegularMatrix<Result, T>& result)
	{
		dispatch_accum_cols(values, J, result, lmat::sum_folder<T>());
	}

	template<typename TI, class JSubs, typename T, class Values, class Result>
	inline void dispatch_max_cols(
			const IEWiseMatrix<Values, T>& values,
			const IEWiseMatrix<JSubs, TI>& J,
			IRegularMatrix<Result, T>& result)
	{
		dispatch_accum_cols(values, J, result, lmat::maximum_folder<T>());
	}

	template<typename TI, class JSubs, typename T, class Values, class Result>
	inline void dispatch_min_cols(
			const IEWiseMatrix<Values, T>& values,
			const IEWiseMatrix<JSubs, TI>& J,
			IRegularMatrix<Result, T>& result)
	{
		dispatch_accum_cols(values, J, result, lmat::minimum_folder<T>());
	}


	template<typename TI, class ISubs, typename T, class Values, class Result, class Folder>
	inline typename std::enable_if<
		lmat::supports_linear_macc<ISubs>::value &&
		is_percol_contiguous<Values>::value &&
		is_percol_contiguous<Result>::value,
	void>::type
	dispatch_accum_rows(
			const IEWiseMatrix<Values, T>& values,
			const IEWiseMatrix<ISubs, TI>& I,
			IRegularMatrix<Result, T>& result,
			const Folder& folder)
	{
		const Values& v = values.derived();
		Result& r = result.derived();

		const index_t m = v.nrows();
		const index_t n = v.ncolumns();
		const index_t K = r.nrows();

		check_arg( I.nelems() == m, "The sizes of J and values are inconsistent" );
		check_arg( r.ncolumns() == n, "The numbers of columns in values and result are inconsistent." );

		auto rd_i = lmat::make_vec_accessor(lmat::atags::scalar(), in_(I.derived()));

		for (index_t j = 0; j < n; ++j)
		{
			const T *vj = v.ptr_col(j);
			T *rj = r.ptr_col(j);

			for (index_t i = 0; i < m; ++i)
			{
				index_t ci = rd_i.scalar(i);
				if (ci >= 0 && ci < K)
				{
					folder.fold(rj[ci], vj[i]);
				}
			}
		}
	}


	template<typename TI, class ISubs, typename T, class Values, class Result>
	inline void dispatch_sum_rows(
			const IEWiseMatrix<Values, T>& values,
			const IEWiseMatrix<ISubs, TI>& I,
			IRegularMatrix<Result, T>& result)
	{
		dispatch_accum_rows(values, I, result, lmat::sum_folder<T>());
	}

	template<typename TI, class ISubs, typename T, class Values, class Result>
	inline void dispatch_max_rows(
			const IEWiseMatrix<Values, T>& values,
			const IEWiseMatrix<ISubs, TI>& I,
			IRegularMatrix<Result, T>& result)
	{
		dispatch_accum_rows(values, I, result, lmat::maximum_folder<T>());
	}

	template<typename TI, class ISubs, typename T, class Values, class Result>
	inline void dispatch_min_rows(
			const IEWiseMatrix<Values, T>& values,
			const IEWiseMatrix<ISubs, TI>& I,
			IRegularMatrix<Result, T>& result)
	{
		dispatch_accum_rows(values, I, result, lmat::minimum_folder<T>());
	}

}

#endif 
