/**
 * @file kmeans_impl.h
 *
 * @brief Internal implementation of K-means
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef DOLPHIN_KMEANS_IMPL_H_
#define DOLPHIN_KMEANS_IMPL_H_

#include <dolphin/common/common_base.h>
#include <dolphin/common/properties.h>

#include <dolphin/common/import_lmat.h>

namespace dolphin { namespace internal {


	template<typename T, typename TI, typename TL,
		class Data, class Centers, class Labels, class Counts, class Monitor>
	void kmeans_impl(
			const Data& data,
			Centers& centers,
			Labels& labels,
			Counts& counts,
			Monitor& monitor,
			const size_t max_iter,
			const T tol)
	{
		const index_t n = data.ncolumns();
		const index_t K = centers.ncolumns();

	}


} }

#endif
