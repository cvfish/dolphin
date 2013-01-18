/**
 * @file kmeans.h
 *
 * @brief K-means algorithm
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef DOLPHIN_KMEANS_H_
#define DOLPHIN_KMEANS_H_

#include "internal/kmeans_impl.h"

namespace dolphin
{

	template<typename T=double>
	class kmeans
	{
	public:
		simple_property<size_t> max_iters;
		simple_property<T> tol;

	public:
		kmeans(size_t K_)
		: max_iters (100,       require_gt_<size_t>(0), "maxiters must be positive")
		, tol       (T(1.0e-6), require_gt_<T>(0),      "tol must be positive")
		{ }

		template<typename TL, typename TI,
			class Data, class Centers, class Labels, class Counts, class Mon>
		void run(const IRegularMatrix<Data, T>& data,
				IRegularMatrix<Centers, T>& centers,
				IRegularMatrix<Labels, TL>& labels,
				IRegularMatrix<Counts, TI>& counts,
				Mon& monitor)
		{
			static_assert(is_percol_contiguous<Data>::value, "data must be percol-contiguous");
			static_assert(is_percol_contiguous<Centers>::value, "centers must be percol-contiguous");
			static_assert(supports_linear_index<Labels>::value, "labels must support linear indexing");
			static_assert(supports_linear_index<Counts>::value, "counts must support linear indexing");

			check_arg(data.nrows() == centers.nrows(),
					"The sample dimensions in data and centers are inconsistent.");

			const index_t n = data.ncolumns();
			const index_t K = centers.ncolumns();

			check_arg(is_vector(labels) && labels.nelems() == n, "The size of labels is invalid.");
			check_arg(is_vector(counts) && counts.nelems() == K, "The size of counts is invalid.");

			internal::kmeans_impl(
					data.derived(),
					centers.derived(),
					labels.derived(),
					counts.derived(),
					monitor, max_iters, tol);
		}
	};


}

#endif
