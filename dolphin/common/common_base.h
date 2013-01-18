/**
 * @file common_base.h
 *
 * @brief Basic definitions to be used over the library
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef DOLPHIN_COMMON_BASE_H_
#define DOLPHIN_COMMON_BASE_H_

#include <light_mat/common/basic_defs.h>

#define DOLPHIN_ENSURE_INLINE LMAT_ENSURE_INLINE

namespace dolphin
{
	// Import of basic types

	using lmat::uint64_t;
	using lmat::int64_t;
	using lmat::uint32_t;
	using lmat::int32_t;
	using lmat::uint16_t;
	using lmat::int16_t;
	using lmat::uint8_t;
	using lmat::uint8_t;

	using lmat::index_t;
	using lmat::size_t;

	using lmat::noncopyable;

	using lmat::invalid_argument;
	using lmat::invalid_operation;

	using lmat::check_arg;
}

#endif
